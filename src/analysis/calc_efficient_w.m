%{

    Compute the efficient weighting matrix needed for running the second-stage efficient GMM.

%}

function [W_matrix] = calc_efficient_w(theta,gmm_function)
    
    % Evaluate model at estimated parameter vector *theta*.
    [~,G_obs_1, G_obs_2, G_obs_3, G_obs_4, G_obs_5, ...
    G_obs_6] = gmm_function(theta);
    
    % Expand micro moment conditions.
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
    
    % Stack variances of moment conditions.
    Omega = blkdiag(Omega_1, Omega_2, Omega_3, Omega_4, Omega_5, Omega_61, Omega_62, Omega_63,Omega_64);
    % Stack in array for easier daignostics in case of rank deficiency.
    Omega_all = {Omega_1,Omega_2,Omega_3,Omega_4,Omega_5,Omega_61,Omega_62,Omega_63,Omega_64};
      
    % Get means of moment matrix.
    g_bar_1 = mean(G_obs_1,1);
    g_bar_2 = mean(G_obs_2,1);
    g_bar_3 = mean(G_obs_3,1);
    g_bar_4 = mean(G_obs_4,1);
    g_bar_5 = mean(G_obs_5,1);
    g_bar_61 = mean(G_obs_61,1);
    g_bar_62 = mean(G_obs_62,1);
    g_bar_63 = mean(G_obs_63,1);
    g_bar_64 = mean(G_obs_64,1);

    % Compute outer product of mean.
    G_bar_1 = g_bar_1' * g_bar_1;
    G_bar_2 = g_bar_2' * g_bar_2;
    G_bar_3 = g_bar_3' * g_bar_3;
    G_bar_4 = g_bar_4' * g_bar_4;
    G_bar_5 = g_bar_5' * g_bar_5;
    G_bar_61 = g_bar_61' * g_bar_61;
    G_bar_62 = g_bar_62' * g_bar_62;
    G_bar_63 = g_bar_63' * g_bar_63;
    G_bar_64 = g_bar_64' * g_bar_64;
    
    % Stack in block diagonal matrix.
    G_bar = blkdiag(G_bar_1, G_bar_2, G_bar_3, G_bar_4, G_bar_5, G_bar_61, G_bar_62, G_bar_63, G_bar_64);
   
    % Construct Omega matrix for xi,nu and zeta-part.
    W = Omega - G_bar;
    % Optimal weighting matrix is inverse of variance of moment conditions.
    W_matrix = inv(W);

    
    % Do rank deficiency check for moment matrices.
    rank_check_Omega = zeros(size(Omega_all,2),1);
    for i=1:size(Omega_all,2)
        Omega_aux = Omega_all{i};
        rank_check_Omega(i) = size(Omega_aux,2) - rank(Omega_aux);
        if rank_check_Omega(i)>0
            fprintf('Rank deficiency in variance matrix block %d: %d dependent columns.\n',i,rank_check_Omega(i));
        end
    end
    % Print diagnostics of efficient weighting matrix.
    rank_check_W_eff = size(W_matrix,2) - rank(W_matrix);
        if rank_check_W_eff>0
            fprintf('Rank deficiency in efficient weighting matrix %d: %d dependent columns.\n',i,rank_check_W_eff);
        end  
end






