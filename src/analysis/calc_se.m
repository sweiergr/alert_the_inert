%{
    Function to calculate standard errors.
    This function takes the step size and the type of perturbation
    (additive or multiplicative) as additional function arguments.
    
%} 

function [varcov, se] = calc_se(theta,gmm_function, WMatrix, ptb_size, ptb_mode)

    %% Set global variables and auxiliary parameters.
    global J T model_type run_estimation SE_calc
    
    % Define model type file suffix.
    SE_calc = 0;
    file_suffix = ['mod',num2str(model_type)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load(project_paths('OUT_ANALYSIS',['thetalinear_',file_suffix,'.mat']),'theta_linear');
    load(project_paths('OUT_ANALYSIS',['expmvalold_',file_suffix,'.mat']),'expmvalold');

    % Define model-type specific file suffix.
    % Extract number of parameters.
    K_l = size(theta_linear,1);
    K_nl = size(theta,1);
    K = K_nl + K_l;
    
    %% Start computation of standard errors.
    % Evaluate model at estimated parameter vector *theta*.
    run_estimation = 1;
    % [~,G_obs_xi, G_obs_zeta, G_obs_mm_sw_ratio, G_obs_mm_pcw_1, G_obs_mm_pcw_2, G_obs_mm, G_obs_mm_adv] = gmm_function(theta);
    [~,G_obs_1, G_obs_2, G_obs_3, G_obs_4, G_obs_5, G_obs_6] = gmm_function(theta);
    run_estimation = 0;
    % Expand micro moment conditions for contract choice with constant and advertising expenditure.
    % Unpack micromoments array.
    G_obs_61 = G_obs_6{1};
    G_obs_62 = G_obs_6{2};
    G_obs_63 = G_obs_6{3};
    G_obs_64 = G_obs_6{4};

    % Get number of observations for different sets of moments.
    n1 = size(G_obs_1,1);
    n2 = size(G_obs_2,1);
    n3 = size(G_obs_3,1);
    n4 = size(G_obs_4,1);
    n5 = size(G_obs_5,1);
    n61 = size(G_obs_61,1);
    n62 = size(G_obs_62,1);
    n63 = size(G_obs_63,1);
    n64 = size(G_obs_64,1);
    n6 = n61+n62+n63+n64;
    % Get number of moments for different sets of moments.
    l1 = size(G_obs_1,2);
    l2 = size(G_obs_2,2);
    l3 = size(G_obs_3,2);
    l4 = size(G_obs_4,2);
    l5 = size(G_obs_5,2);
    l61 = size(G_obs_61,2);
    l62 = size(G_obs_62,2);
    l63 = size(G_obs_63,2);
    l64 = size(G_obs_64,2);
    l6 = l61+l62+l63+l64;

    % Select correct number of observations here.
    n_obs = n6;
    % Select number of moment conditions. Not much wiggle room here.
    n_mc_total = l1 + l2 + l3 + l4 + l5 + l6; 
    
    % Construct Omega matrices for the six different sets of moments.
    Omega_1 = (1/n1) .* (G_obs_1' * G_obs_1);
    Omega_2 = (1/n2) .* (G_obs_2' * G_obs_2);
    Omega_3 = (1/n3) .* (G_obs_3' * G_obs_3);
    Omega_4 = (1/n4) .* (G_obs_4' * G_obs_4);
    Omega_5 = (1/n5) .* (G_obs_5' * G_obs_5);
    Omega_61 = (1/n61) .* (G_obs_61' * G_obs_61);
    Omega_62 = (1/n62) .* (G_obs_62' * G_obs_62);
    Omega_63 = (1/n63) .* (G_obs_63' * G_obs_63);
    Omega_64 = (1/n64) .* (G_obs_64' * G_obs_64);


    % Stack outer products of moment conditions.
    Omega = blkdiag(Omega_1, Omega_2, Omega_3, ...
                    Omega_4, Omega_5, ...
                    Omega_61, Omega_62, Omega_63,Omega_64);
    
    % Initialize matrices for derivatives parts.
    % Rows: number of moment conditions, columns: number of parameters.
    Q_1 = zeros(l1,K);
    Q_2 = zeros(l2,K);
    Q_3 = zeros(l3,K);
    Q_4 = zeros(l4,K);
    Q_5 = zeros(l5,K);
    Q_61 = zeros(l61,K);
    Q_62 = zeros(l62,K);
    Q_63 = zeros(l63,K);
    Q_64 = zeros(l64,K);
    
    % Define index threshold for advertising parameters.
    adv_idx = length(theta) + length(theta_linear) -2;
    % Initiate actual SE computation.
    SE_calc = 1;
    % Loop over parameters.
    for ii=1:K
        % Reload linear parameter vector.
        load(project_paths('OUT_ANALYSIS',['thetalinear_',file_suffix,'.mat']),'theta_linear');
        % Because advertising moments are "less continuous" we use a slightly larger perturbation for these parameters.
        if ii>adv_idx
            ptb_size = 1E-2;
        end
        
        % Perturb "linear" parameters.
        if ii<=K_l
            % Perturb parameter ii.
            hh = ptb_size; % * theta_linear(ii,1);
            h = (1:K_l==ii)';
            if isequal(ptb_mode,'add')
                theta_linear = theta_linear - hh .* h;
            elseif isequal(ptb_mode,'mult')
                theta_linear = theta_linear - hh .* (h .* theta_linear);
            end

            theta_grad_b = theta;
            theta_grad_f = theta; 
        % Perturb "nonlinear" parameter
        elseif ii>K_l
            hh = ptb_size; % * theta(ii,1);
            h = (1:K_nl==ii-K_l)';
            if isequal(ptb_mode,'add')
                theta_grad_b = theta - hh .* h;
                theta_grad_f = theta + hh .* h;
            elseif isequal(ptb_mode,'mult')
                theta_grad_b = theta - hh .* (h .* theta);
                theta_grad_f = theta + hh .* (h .* theta);
            end
        end
        % Save perturbed (if necessary) linear parameter vector to file.
        save(project_paths('OUT_ANALYSIS',['thetalinear_SE_aux_',file_suffix,'.mat']),'theta_linear');

        % Evaluate model at perturbed parameter vector.
        [~, G_1_ptb, G_2_ptb, G_3_ptb, G_4_ptb, G_5_ptb,  G_6_ptb] = gmm_function(theta_grad_b);
        
        G_61_ptb = G_6_ptb{1};
        G_62_ptb = G_6_ptb{2};
        G_63_ptb = G_6_ptb{3};
        G_64_ptb = G_6_ptb{4};
       
        load(project_paths('OUT_ANALYSIS',['thetalinear_',file_suffix,'.mat']),'theta_linear');
        if ii<=K_l
            if isequal(ptb_mode,'add')
                theta_linear = theta_linear + hh .* h;
            elseif isequal(ptb_mode,'mult')
                theta_linear = theta_linear + hh .* (h .* theta_linear);
            end
       end
       % Save linear parameter vector to aux file.
       save(project_paths('OUT_ANALYSIS',['thetalinear_SE_aux_',file_suffix,'.mat']),'theta_linear');
       % Call model file with perturbed parameters. 
       [~, G_1_ptf, G_2_ptf, G_3_ptf, G_4_ptf, G_5_ptf, G_6_ptf] = gmm_function(theta_grad_f);
        G_61_ptf = G_6_ptf{1};
        G_62_ptf = G_6_ptf{2};
        G_63_ptf = G_6_ptf{3};
        G_64_ptf = G_6_ptf{4};
        
        % Get derivative approximation using finite-differences.
        G_1_prime = (2*hh)^-1 .* (G_1_ptf - G_1_ptb);
        G_2_prime = (2*hh)^-1 .* (G_2_ptf - G_2_ptb);
        G_3_prime = (2*hh)^-1 .* (G_3_ptf - G_3_ptb);
        G_4_prime = (2*hh)^-1 .* (G_4_ptf - G_4_ptb);
        G_5_prime = (2*hh)^-1 .* (G_5_ptf - G_5_ptb);
        G_61_prime = (2*hh)^-1 .* (G_61_ptf - G_61_ptb);
        G_62_prime = (2*hh)^-1 .* (G_62_ptf - G_62_ptb);
        G_63_prime = (2*hh)^-1 .* (G_63_ptf - G_63_ptb);
        G_64_prime = (2*hh)^-1 .* (G_64_ptf - G_64_ptb);
       
        % Compute expected value of derivatives.
        G_1_prime_mean = mean(G_1_prime,1); 
        G_2_prime_mean = mean(G_2_prime,1); 
        G_3_prime_mean = mean(G_3_prime,1);
        G_4_prime_mean = mean(G_4_prime,1);
        G_5_prime_mean = mean(G_5_prime,1);
        G_61_prime_mean = mean(G_61_prime,1);
        G_62_prime_mean = mean(G_62_prime,1);
        G_63_prime_mean = mean(G_63_prime,1);
        G_64_prime_mean = mean(G_64_prime,1);

        Q_1(:,ii) = G_1_prime_mean;
        Q_2(:,ii) = G_2_prime_mean;
        Q_3(:,ii) = G_3_prime_mean;
        Q_4(:,ii) = G_4_prime_mean;
        Q_5(:,ii) = G_5_prime_mean;
        Q_61(:,ii) = G_61_prime_mean;
        Q_62(:,ii) = G_62_prime_mean;
        Q_63(:,ii) = G_63_prime_mean;
        Q_64(:,ii) = G_64_prime_mean;
    end
    % Stack derivatives for all sets of moment conditions.
    Q = [Q_1; Q_2; Q_3; Q_4; Q_5; ...
         Q_61; Q_62; Q_63; Q_64 ...
         ];
    
    %% Put formula for variance-covariance matrix together.
    varcov = ((Q'*WMatrix*Q) \ (Q'*WMatrix*Omega*WMatrix*Q)) / (Q'*WMatrix*Q) ./ n_obs;
    se = sqrt(diag(varcov));
    % Reset flag for SE computation to zero.
    SE_calc = 0;
end
