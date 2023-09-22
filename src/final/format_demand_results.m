%{
    Write demand estimates for preferenecs to LaTeX table.	
    This file now handles all model specifications.

%}


% Define list of models to format.
% It seems like models 1,4,5,6 are most promising at this point. 
    % Not sure what 2 and 3 add.
    % 1: both psi and kappa are homogeneous.
    % These models probably don't add much for now.
    % 2: psi is homogeneous, kappa is heterogeneous.
    % 3: psi is heterogeneous, kappa is homogeneous.
    % 4: both psi and kappa are heterogeneous.
    % 5: only search frictions, i.e., no switching costs
    % 6: only switching costs (capturing everything in reduced form), i.e., no search costs from using PCW.
model_list_format = [1,4,5,6];
% Define whether to two efficient two-step GMM or just first-step.
% rows indicate model to esitmate.
% first column indicates first step GMM, second column efficient two-step GMM.
gmm_list_format = [1,2];
for model_idx_format=1:length(model_list_format)
    for gmm_idx_format=1:length(gmm_list_format)
        clearvars -except model_list_format gmm_list_format model_idx_format gmm_idx_format
        % Set which model type to run.
        model_type = model_list_format(model_idx_format);
        % Define model type file suffix.
        file_suffix = [num2str(gmm_idx_format)','_','mod',num2str(model_type)];
        fprintf('Formatting results for model specification %d and GMM stage %d...\n', model_type,gmm_idx_format);
        % Load workspace from GMM estimation.
        load(project_paths('OUT_ANALYSIS',['postestimation_workspace_', file_suffix,'.mat']));
        file_suffix = [num2str(gmm_idx_format)','_','mod',num2str(model_type)];
        % Define Row labels.
        Varnames={'Price'; 'Incumbent'; 'Green Electricity'; 'Switching Costs'}; 
        if gmm_idx_format==2
            theta_opt = theta_opt_2;
            SE = SE_2;
            p_values = p_values_2;
        end
        % Number of nonlinear parameters.
        k_lin = length(theta_linear);
        k_nl = length(theta_opt);
        %% SW: THIS PART DEPENDS ON MODEL SPECIFICATION!
        % Arrange matrix of coefficients and let 0 denote parameters that are not
        % estimated
        if model_type==1
            coeff_mean = [theta_opt(4); theta_linear(1); theta_linear(2); theta_opt(8)];
            coeff_sigma = [zeros(2,1); theta_opt(5); zeros(1,1)];
            coeff_age = [zeros(1,1); theta_opt(6); zeros(2,1)];
            coeff_income = [theta_opt(7); zeros(3,1)];
            coeff_matrix = [coeff_mean, coeff_sigma, coeff_age, coeff_income];
            % matrix of standard errors.
            se_mean = [SE(k_lin+4); SE(1); SE(2); SE(k_lin+8)];
            se_sigma = [zeros(2,1); SE(k_lin+5); zeros(1,1)];
            se_age = [zeros(1,1); SE(k_lin+6); zeros(2,1)];
            se_income = [SE(k_lin+7); zeros(3,1)];
            se_matrix = [se_mean, se_sigma, se_age, se_income];

            % Matrix of p-values.
            pval_mean = [p_values(k_lin+4); p_values(1); p_values(2);p_values(k_lin+8)];
            pval_sigma = [NaN(2,1); p_values(k_lin+5); NaN(1,1)];
            pval_age = [NaN(1,1); p_values(k_lin+6); NaN(2,1)];
            pval_income = [p_values(k_lin+7); NaN(3,1)];
            pval_matrix = [pval_mean, pval_sigma, pval_age, pval_income];
            % Define vector to indicate which parameter interactiosn are estimated.
            var_check = [1001;1010;1100;1000];
        elseif model_type==4 % large model with both heterogeneous switching costs and PCW search costs
            coeff_mean = [theta_opt(5); theta_linear(1); theta_linear(2); theta_opt(9)];
            coeff_sigma = [zeros(2,1); theta_opt(6); zeros(1,1)];
            coeff_age = [zeros(1,1); theta_opt(7); zeros(1,1);theta_opt(10)];
            coeff_income = [theta_opt(8); zeros(3,1)];
            coeff_matrix = [coeff_mean, coeff_sigma, coeff_age, coeff_income];
            % matrix of standard errors.
            se_mean = [SE(k_lin+5); SE(1); SE(2); SE(k_lin+9)];
            se_sigma = [zeros(2,1); SE(k_lin+6); zeros(1,1)];
            se_age = [zeros(1,1); SE(k_lin+7); zeros(1,1);SE(k_lin+10)];
            se_income = [SE(k_lin+8); zeros(3,1)];
            se_matrix = [se_mean, se_sigma, se_age, se_income];
            % Matrix of p-values.
            pval_mean = [p_values(k_lin+5); p_values(1); p_values(2);p_values(k_lin+9)];
            pval_sigma = [NaN(2,1); p_values(k_lin+6); NaN(1,1)];
            pval_age = [NaN(1,1); p_values(k_lin+7); NaN(1,1);p_values(k_lin+10)];
            pval_income = [p_values(k_lin+8); NaN(3,1)];
            pval_matrix = [pval_mean, pval_sigma, pval_age, pval_income];
            % Define vector to indicate which parameter interactiosn are estimated.
            var_check = [1001;1010;1100;1010];
      elseif model_type==5 % restricted model without switching costs
            coeff_mean = [theta_opt(4); theta_linear(1); theta_linear(2); zeros(1,1)];
            coeff_sigma = [zeros(2,1); theta_opt(5); zeros(1,1)];
            coeff_age = [zeros(1,1); theta_opt(6); zeros(2,1)];
            coeff_income = [theta_opt(7); zeros(3,1)];
            coeff_matrix = [coeff_mean, coeff_sigma, coeff_age, coeff_income];
            % matrix of standard errors.
            se_mean = [SE(k_lin+4); SE(1); SE(2); zeros(1,1)];
            se_sigma = [zeros(2,1); SE(k_lin+5); zeros(1,1)];
            se_age = [zeros(1,1); SE(k_lin+6); zeros(2,1)];
            se_income = [SE(k_lin+7); zeros(3,1)];
            se_matrix = [se_mean, se_sigma, se_age, se_income];
            % Matrix of p-values.
            pval_mean = [p_values(k_lin+4); p_values(1); p_values(2);NaN(1,1)];
            pval_sigma = [NaN(2,1); p_values(k_lin+5); NaN(1,1)];
            pval_age = [NaN(1,1); p_values(k_lin+6); NaN(2,1)];
            pval_income = [p_values(k_lin+7); NaN(3,1)];
            pval_matrix = [pval_mean, pval_sigma, pval_age, pval_income];
            % Define vector to indicate which parameter interactiosn are estimated.
            var_check = [1001;1010;1100;1000];
      elseif model_type==6 % restricted model without PCW search costs
            coeff_mean = [theta_opt(1); theta_linear(1); theta_linear(2); theta_opt(5)];
            coeff_sigma = [zeros(2,1); theta_opt(2); zeros(1,1)];
            coeff_age = [zeros(1,1); theta_opt(6); zeros(2,1)];
            coeff_income = [theta_opt(4); zeros(3,1)];
            coeff_matrix = [coeff_mean, coeff_sigma, coeff_age, coeff_income];
            % matrix of standard errors.
            se_mean = [SE(k_lin+1); SE(1); SE(2); SE(k_lin+5)];
            se_sigma = [zeros(2,1); SE(k_lin+2); zeros(1,1)];
            se_age = [zeros(1,1); SE(k_lin+3); zeros(2,1)];
            se_income = [SE(k_lin+4); zeros(3,1)];
            se_matrix = [se_mean, se_sigma, se_age, se_income];
            % Matrix of p-values.
            pval_mean = [p_values(k_lin+1); p_values(1); p_values(2);p_values(k_lin+5)];
            pval_sigma = [NaN(2,1); p_values(k_lin+2); NaN(1,1)];
            pval_age = [NaN(1,1); p_values(k_lin+3); NaN(2,1)];
            pval_income = [p_values(k_lin+4); NaN(3,1)];
            pval_matrix = [pval_mean, pval_sigma, pval_age, pval_income];
            % Define vector to indicate which parameter interactiosn are estimated.
            var_check = [1001;1010;1100;1000];
      end
        sig_stars_legend = {'','*','**','***'};
        sig_stars_matrix = 1 .* (pval_matrix>0.1 | isnan(pval_matrix)) + 2 .* (pval_matrix<=0.1 & pval_matrix>0.05) + 3 .* (pval_matrix<=0.05 & pval_matrix>0.01) +  4 .* (pval_matrix<0.01);
        sig_stars = sig_stars_legend(sig_stars_matrix);

        % Open file handle to write table to..
        fid = fopen(project_paths('OUT_TABLES',['estresultspreferences_',file_suffix,'.tex']),'w');

        fprintf(fid,'\\begin{tabular}{ lcccc  } \\toprule \n');
        %fprintf(fid, '& Mean & St.Dev. $\\sigma$ & Children & College & Income  \\\\ \n');
        fprintf(fid, '& Mean & Sigma & Senior &  Income  \\\\ \n');
        fprintf(fid,'\\midrule \n');
        % Loop over all variables.
        for k=1:length(Varnames)  
            % Row for linear coefficient only.
            if var_check(k,1)==1000
                fprintf(fid, '%s & $%8.4f$%s & &  &   \\\\ \n ', Varnames{k,1}, coeff_matrix(k,1), sig_stars{k,1} );   
                fprintf(fid, '		& $(%8.4f)$ & &  &   \\\\ \n', se_matrix(k,1));

            % Row for mean and senior interaction.
            elseif var_check(k,1)==1010
                fprintf(fid, '%s & $%8.4f$%s & & $%8.4f$%s &     \\\\ \n ', Varnames{k,1}, coeff_matrix(k,1), sig_stars{k,1}, coeff_matrix(k,3), sig_stars{k,3});   
                fprintf(fid, '	 & $(%8.4f)$ & & $(%8.4f)$ &  \\\\ \n', se_matrix(k,1), se_matrix(k,3));

            % Row for mean and sigma (random coefficient).
            elseif var_check(k,1)==1100
                fprintf(fid, '%s & $%8.4f$%s & $%8.4f$%s &  & \\\\ \n ', Varnames{k,1}, coeff_matrix(k,1), sig_stars{k,1}, coeff_matrix(k,2), sig_stars{k,2});   
                fprintf(fid, '   & $(%8.4f)$ & $(%8.4f)$ &  & \\\\ \n', se_matrix(k,1), se_matrix(k,2));

            % Row for mean and income interaction.    
            elseif var_check(k,1)==1001
                fprintf(fid, '%s & $%8.4f$%s&  & & $%8.4f$%s \\\\ \n ', Varnames{k,1}, coeff_matrix(k,1), sig_stars{k,1}, coeff_matrix(k,4), sig_stars{k,4});   
                fprintf(fid, '   & $(%8.4f)$&  & & $(%8.4f)$ \\\\ \n', se_matrix(k,1), se_matrix(k,4));
           end
            % Write bottom rule
            if k==length(Varnames)
                fprintf(fid, '\\midrule \n ');
            end
        end
        if gmm_idx_format==1
            fprintf(fid,'\\multicolumn{5}{p{12.5cm}}{\\footnotesize{\\textit{Notes: Results for preference parameters from estimating the demand model using GMM with block-diagonal 2SLS weighting matrix. Standard errors in parentheses. *,**,*** denote significance at the 10, 5 and 1 percent level respectively.}}}\\\\ \n');
        elseif gmm_idx_format==2
            fprintf(fid,'\\multicolumn{5}{p{12.5cm}}{\\footnotesize{\\textit{Notes: Results for preference parameters from estimating the demand model using efficient two-step GMM. Standard errors in parentheses. *,**,*** denote significance at the 10, 5 and 1 percent level respectively.}}}\\\\ \n');
        end
            %fprintf(fid,'\\multicolumn{6}{l}{\\footnotesize{\\textit{}}}\\\\ \n')
        fprintf(fid, '\\bottomrule \n');
        % Close tabular part.
        fprintf(fid,  '\\end{tabular}\n');
        % Close file handle.
        fclose(fid);
    end % end loop over GMM stages.   
end % end loop over model specifications.