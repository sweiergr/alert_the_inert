%{ 

Load post-estimation workspaces and run counterfactuals for each of the models discussed in the main text.
This file generates the results for the counterfactual table (Table 3 in the main text).
Formatting of the results is done in *src/final/format_counterfactuals.py*.

%}

%% General preamble clearing memory and creating log-file.
clc
clear
mex('-R2018a',project_paths('IN_ANALYSIS','contrMap_ES.c'));

%% Set global variables.
global T J JF NS_cons NS_p_eps ...
    run_counterfactuals SE_calc run_estimation ...
    cf_subsidy_switch subsidy_switch ...
    cf_perfect_info kappa_cut ...
    cf_reg_campaign cf_inc_only_all debug_model model_type %g_bar_all theta_all model_type
%% Load relevant workspace: depends on which spec we want to use for counterfactuals.
debug_model = 0;
% Indicator for whether switching costs work on firm-level (0) or
% contract-level (1).
% Contract-level requires resolving of model for counterfactuals, where
% several contracts are offered.
contract_sc = 1;

% This indicates the additional surcharge for a kwh of green electricity in
% EUR. This is only used for calculating indifference markup for some of the
% counterfactuals where green electricity contracts are available.
% The results are not affected significantly by choosing this within a reasonable range.
green_surcharge_factor = 0.005;
% Define list of models to use for counterfactual.
    % Not sure what 2 and 3 add.
    % 1 (baseline model): both psi and kappa are homogeneous.
        % These models probably don't add much for now.
        % 2: psi is homogeneous, kappa is heterogeneous.
        % 3: psi is heterogeneous, kappa is homogeneous.
    % 4 (model extension reported in main text): both psi and kappa are heterogeneous.
    % 5 (restricted model 1): only search frictions, i.e., no switching costs
    % 6 (restricted model 2): only switching costs (capturing everything in reduced form), i.e., no search costs from using PCW.
model_list_cf = [1,4,5,6];
% Define whether to use first-stage of efficient second-stage results for
% each model type.
% first column indicates first step GMM, second column efficient two-step GMM.
gmm_type_cf = [1; ... for baseline model 1
            1; ... for bigger model 2 (currently not used)
            1; ... for bigger model 3 (currently not used)
            1; ... for bigger model with both frictions heterogeneous across age groups
            0; ... for model without switching cost (strawman number 1)
            0]; %for model without PCW search costs (strawman number 2)
% Loop over different models.
for model_idx=1:length(model_list_cf)
%% Set which model type to run.
model_type = model_list_cf(model_idx);
% Define model type file suffix.
file_suffix = ['mod',num2str(model_type)];
% Whether to use first stage of second stage efficient estimates.
efficient_GMM = gmm_type_cf(model_type);
if efficient_GMM==0
    file_suffix_full = [num2str(efficient_GMM+1),'_','mod',num2str(model_type)];
elseif efficient_GMM==1
    file_suffix_full = [num2str(efficient_GMM+1),'_','mod',num2str(model_type)];
end
% Load relevant workspace and mean utilities and linear parameters.
if efficient_GMM==0
    load(project_paths('OUT_ANALYSIS',['postestimation_workspace_1_',file_suffix]));
    % Load linear parameters from separate file.
    %theta_linear_backup = theta_linear;
    load(project_paths('OUT_ANALYSIS',['theta_lin_opt_1_',file_suffix]));
    theta_linear = theta_lin_opt_1;
    save(project_paths('OUT_ANALYSIS',['theta_linear_',file_suffix,'.mat']),'theta_linear');
    % Load relevant deltas.
    load(project_paths('OUT_ANALYSIS',['expmvalold_opt_1_',file_suffix]));
    expmvalold = expmval_opt_1;
    save(project_paths('OUT_ANALYSIS',['expmvalold_',file_suffix,'.mat']),'expmvalold');
    % This is used as a basis for the incumbent only selling all contracts
    % below.
    expmvalold_backup = expmvalold;
    expmvalold_track = expmvalold_backup;
    theta_linear_backup = theta_linear;
    % Update anoynmous model function and desired convergence criteria (if relevant).
    % This is safety check #1: If we run the contraction mapping nothing
    % should change anymore.
    % Specify convergence criteria for contraction mapping routine.
    pred_shares_track = market_shares;
    run_estimation=1;
    tol_delta_c = 1E-6;
    max_iter_c = 1000;
    gmm_function = @(theta) model_myopic(theta, cdid, cdindex, market_shares, s_init_long, avg_churn_rates, mm_data_reshaped, pcw_use_micro, pcw_inst_micro, pcw_use_macro, adv_data, adv_data_mm, data_demo, X_reg, Z_matrix, draws_p_normalized, draws_eps,  preference_shocks,awareness_shocks,WMatrix, tol_delta_c,max_iter_c);
    % Do another model function call at optimal parameter values to ensure convergence
    % and get predicted market shares, churn rates, and PCW usage.
    [~,~,~,~,~,~,~, ...
        pred_shares1d, pred_churnAgg1d, ...
        cSurplusPCW_obs,cSurplusNoPCW_obs, ...
        pred_PCWAgg1d, ...
        mu_obs,price_coeff_obs, ...
        pred_CCPAvg,~,~,cSurplusI_obs] = gmm_function(theta_opt);
    % Sanity check code.
    pred_shares_track = [pred_shares_track pred_shares1d];
    % Check that linear parameters and mean utilities have not changed.
    load(project_paths('OUT_ANALYSIS',['thetalinear_',file_suffix,'.mat']),'theta_linear');
    load(project_paths('OUT_ANALYSIS',['expmvalold_',file_suffix,'.mat']),'expmvalold');
    save(project_paths('OUT_ANALYSIS',['theta_linear_',file_suffix,'.mat']),'theta_linear');
    save(project_paths('OUT_ANALYSIS',['expmvalold_',file_suffix,'.mat']),'expmvalold');
    conv_check_meanutil_1 = norm(expmval_opt_1-expmvalold);
    conv_check_meanutil_long = expmval_opt_1-expmvalold;
    conv_check_idx = (abs(conv_check_meanutil_long)>0.1);
    conv_check_delta_wide = reshape(conv_check_meanutil_long,J,T)';
    conv_check_idx_wide = reshape(conv_check_idx,J,T)';
    conv_check_thetalin_1 = norm(theta_lin_opt_1-theta_linear);
    conv_check_pshares = norm(pred_shares_track(12:end,end) - pred_shares_track(12:end,end-1));
    conv_check_pshares_full = norm(pred_shares_track(1:end,end) - pred_shares_track(1:end,end-1));
    fprintf('Safety iteration 1: Norm difference between exponentiated mean utilities obtained after first stage GMM is:\n %6.4f\n',conv_check_meanutil_1);
    fprintf('Safety iteration 1: Norm difference between linear parameters obtained after first stage GMM is:\n %6.4f\n',conv_check_thetalin_1);
    fprintf('Safety iteration 1 (first period excluded): Norm difference between predicted shares obtained after first stage GMM is:\n %6.4f\n',conv_check_pshares);    
    fprintf('Safety iteration 1 (first period included): Norm difference between predicted shares obtained after first stage GMM is:\n %6.4f\n',conv_check_pshares_full);    

    % Create month index.
    month_idx = kron(linspace(1,T,T)',ones(J,1));
    % Track evolution of mean utilities when model function is called
    % several times.
    expmvalold_track = [expmvalold_track, expmvalold];
    expmval_rel_change = expmvalold_track(:,end) ./ expmvalold_track(:,1);
    theta_linear_track = [theta_linear_backup, theta_linear];
    shares_pred_wide = reshape(pred_shares1d,J,T)';
    shares_obs_wide = reshape(market_shares,J,T)';
    % Safety check # 2: Just resolve model and see that nothing changes.
    run_estimation = 0;
    gmm_function = @(theta) model_myopic(theta, cdid, cdindex, market_shares, s_init_long, avg_churn_rates, mm_data_reshaped, pcw_use_micro, pcw_inst_micro, pcw_use_macro, adv_data, adv_data_mm, data_demo, X_reg, Z_matrix, draws_p_normalized, draws_eps,  preference_shocks,awareness_shocks,WMatrix, tol_delta_c,max_iter_c);
    % Do another model function call at optimal parameter values to ensure convergence
    % and get predicted market shares, churn rates, and PCW usage.
    [~,~,~,~,~,~,~, ...
        pred_shares1d, pred_churnAgg1d, ...
        cSurplusPCW_obs,cSurplusNoPCW_obs, ...
        pred_PCWAgg1d, ...
        mu_obs,price_coeff_obs, ...
        pred_CCPAvg,~,~,cSurplusI_obs] = gmm_function(theta_opt);
    % Track evolution of mean utilities when model function is caleed
    % several times.
    load(project_paths('OUT_ANALYSIS',['thetalinear_',file_suffix,'.mat']),'theta_linear');
    load(project_paths('OUT_ANALYSIS',['expmvalold_',file_suffix,'.mat']),'expmvalold');
    expmvalold_track = [expmvalold_track, expmvalold];
    theta_linear_track = [theta_linear_track, theta_linear];
    pred_shares_track = [pred_shares_track pred_shares1d];
    conv_check_meanutil_2 = norm(expmvalold_track(:,end)-expmvalold_track(:,end-1));
    conv_check_thetalin_2 = norm(theta_linear_track(:,end)-theta_linear_track(:,end-1));
    conv_check_pshares_2 = norm(pred_shares_track(1:end,end)-pred_shares_track(1:end,end-1));
    fprintf('Safety iteration 2: Norm difference between exponentiated mean utilities obtained after first stage GMM is:\n %6.4f\n',conv_check_meanutil_2);
    fprintf('Safety iteration 2: Norm difference between linear parameters obtained after first stage GMM is:\n %6.4f\n',conv_check_thetalin_2);
    fprintf('Safety iteration 2 (first period included): Norm difference between predicted shares obtained after first stage GMM is:\n %6.4f\n',conv_check_pshares_2);    
    fprintf('\nEND OF PRE-COUNTERFACTUAL SAFETY CHECKS\n\n');
    % Write optimal parameter vector to new container (to make code general
    % across first and second stage).
    theta_opt_cf = theta_opt;
elseif efficient_GMM==1
    load(project_paths('OUT_ANALYSIS',['postestimation_workspace_2_',file_suffix]));
    % Backup and load linear parameters from separate file.
    %theta_linear_backup = theta_linear;
    load(project_paths('OUT_ANALYSIS',['theta_lin_opt_2_',file_suffix]));
    theta_linear = theta_lin_opt_2;
    save(project_paths('OUT_ANALYSIS',['theta_linear_',file_suffix,'.mat']),'theta_linear');
    % Load relevant deltas.
    load(project_paths('OUT_ANALYSIS',['expmvalold_opt_2_',file_suffix]));
    expmvalold = expmval_opt_2;
    save(project_paths('OUT_ANALYSIS',['expmvalold_',file_suffix,'.mat']),'expmvalold');
    % This is used as a basis for the incumbent only selling all contracts
    % below.
    expmvalold_backup = expmvalold;
    expmvalold_track = expmvalold_backup;
    theta_linear_backup = theta_linear;
    % Update anoynmous model function and desired convergence criteria (if relevant).
    % This is safety check #1: If we run the contraction mapping nothing
    % should change anymore.
    % Specify convergence criteria for contraction mapping routine.
    pred_shares_track = market_shares;
    run_estimation=1;
    tol_delta_c = 1E-6;
    max_iter_c = 1000;
    gmm_function = @(theta) model_myopic(theta, cdid, cdindex, market_shares, s_init_long, avg_churn_rates, mm_data_reshaped, pcw_use_micro, pcw_inst_micro, pcw_use_macro, adv_data, adv_data_mm, data_demo, X_reg, Z_matrix, draws_p_normalized, draws_eps,  preference_shocks,awareness_shocks,WMatrix, tol_delta_c,max_iter_c);
    % Do another model function call at optimal parameter values to ensure convergence
    % and get predicted market shares, churn rates, and PCW usage.
    [~,~,~,~,~,~,~, ...
        pred_shares1d, pred_churnAgg1d, ...
        cSurplusPCW_obs,cSurplusNoPCW_obs, ...
        pred_PCWAgg1d, ...
        mu_obs,price_coeff_obs, ...
        pred_CCPAvg,~,~,cSurplusI_obs] = gmm_function(theta_opt_2);
    % Sanity check code.
    pred_shares_track = [pred_shares_track pred_shares1d];
    
    % Check that linear parameters and mean utilities have not changed.
    load(project_paths('OUT_ANALYSIS',['thetalinear_',file_suffix,'.mat']),'theta_linear');
    load(project_paths('OUT_ANALYSIS',['expmvalold_',file_suffix,'.mat']),'expmvalold');
    save(project_paths('OUT_ANALYSIS',['theta_linear_',file_suffix,'.mat']),'theta_linear');
    save(project_paths('OUT_ANALYSIS',['expmvalold_',file_suffix,'.mat']),'expmvalold');
    conv_check_meanutil_1 = norm(expmval_opt_2-expmvalold);
    conv_check_meanutil_long = expmval_opt_2-expmvalold;
    conv_check_idx = (abs(conv_check_meanutil_long)>0.1);
    conv_check_delta_wide = reshape(conv_check_meanutil_long,J,T)';
    conv_check_idx_wide = reshape(conv_check_idx,J,T)';
    conv_check_thetalin_1 = norm(theta_lin_opt_1-theta_linear);
    conv_check_pshares = norm(pred_shares_track(12:end,end) - pred_shares_track(12:end,end-1));
    conv_check_pshares_full = norm(pred_shares_track(1:end,end) - pred_shares_track(1:end,end-1));
    fprintf('Safety iteration 1: Norm difference between exponentiated mean utilities obtained after second stage GMM is:\n %6.4f\n',conv_check_meanutil_1);
    fprintf('Safety iteration 1: Norm difference between linear parameters obtained after second stage GMM is:\n %6.4f\n',conv_check_thetalin_1);
    fprintf('Safety iteration 1 (first period excluded): Norm difference between predicted shares obtained after second stage GMM is:\n %6.4f\n',conv_check_pshares);    
    fprintf('Safety iteration 1 (first period included): Norm difference between predicted shares obtained after second stage GMM is:\n %6.4f\n',conv_check_pshares_full);    

    % Create month index.
    month_idx = kron(linspace(1,T,T)',ones(J,1));
    % Track evolution of mean utilities when model function is caleed
    % several times.
    expmvalold_track = [expmvalold_track, expmvalold];
    expmval_rel_change = expmvalold_track(:,end) ./ expmvalold_track(:,1);
    theta_linear_track = [theta_linear_backup, theta_linear];
    shares_pred_wide = reshape(pred_shares1d,J,T)';
    shares_obs_wide = reshape(market_shares,J,T)';
    % Safety check # 2: Just resolve model and see that nothing changes.
    run_estimation = 0;
    gmm_function = @(theta) model_myopic(theta, cdid, cdindex, market_shares, s_init_long, avg_churn_rates, mm_data_reshaped, pcw_use_micro, pcw_inst_micro, pcw_use_macro, adv_data, adv_data_mm, data_demo, X_reg, Z_matrix, draws_p_normalized, draws_eps,  preference_shocks,awareness_shocks,WMatrix, tol_delta_c,max_iter_c);
    % Do another model function call at optimal parameter values to ensure convergence
    % and get predicted market shares, churn rates, and PCW usage.
    [~,~,~,~,~,~,~, ...
        pred_shares1d, pred_churnAgg1d, ...
        cSurplusPCW_obs,cSurplusNoPCW_obs, ...
        pred_PCWAgg1d, ...
        mu_obs,price_coeff_obs, ...
        pred_CCPAvg,~,~,cSurplusI_obs] = gmm_function(theta_opt_2);
    % Track evolution of mean utilities when model function is caleed
    % several times.
    load(project_paths('OUT_ANALYSIS',['thetalinear_',file_suffix,'.mat']),'theta_linear');
    load(project_paths('OUT_ANALYSIS',['expmvalold_',file_suffix,'.mat']),'expmvalold');
    expmvalold_track = [expmvalold_track, expmvalold];
    theta_linear_track = [theta_linear_track, theta_linear];
    pred_shares_track = [pred_shares_track pred_shares1d];
    conv_check_meanutil_2 = norm(expmvalold_track(:,end)-expmvalold_track(:,end-1));
    conv_check_thetalin_2 = norm(theta_linear_track(:,end)-theta_linear_track(:,end-1));
    conv_check_pshares_2 = norm(pred_shares_track(1:end,end)-pred_shares_track(1:end,end-1));
    fprintf('Safety iteration 2: Norm difference between exponentiated mean utilities obtained after second stage GMM is:\n %6.4f\n',conv_check_meanutil_2);
    fprintf('Safety iteration 2: Norm difference between linear parameters obtained after second stage GMM is:\n %6.4f\n',conv_check_thetalin_2);
    fprintf('Safety iteration 2 (first period included): Norm difference between predicted shares obtained after second stage GMM is:\n %6.4f\n',conv_check_pshares_2);    
    fprintf('\nEND OF PRE-COUNTERFACTUAL SAFETY CHECKS\n\n');
    % Write optimal parameter vector to new container (to make code general
    % across first and second stage).
    theta_opt_cf = theta_opt_2;
end

%% Set global parameters to manage counterfactuals.
debug_model=0;
run_counterfactuals = 1;
run_estimation==0;
% These are important selectors on how model predictions are computed. Don't delete!
cf_subsidy_switch = 0;
cf_perfect_info = 0;
cf_reg_campaign = 0;

date = datestr(now);
print_details = 1;
format 'short'
% Construct indicator for grid based on OS used.
% Determine whether code is run on local Mac/Windows or supercomputer grid.
os_ind = computer;
% Indicate whether estimation is run on IU grid/Linux computer.
if os_ind(1:3) == 'XXX'
    grid=1;
else
    grid=0;
end
%% Set path to diary log.
if grid==1
    %addpath('functions');
    input_path = 'input_data/';
    output_path = 'output_estimation/';
    output_path_log = output_path;
else
    % Add necessary paths to source libraries.
    %addpath('../library/matlab')
    input_path = project_paths('OUT_DATA');
    output_path = project_paths('OUT_TABLES');
    output_path_log = project_paths('OUT_ANALYSIS');
end
% Write log file to diary.
% Delete existing file if it already exists (to avoid appending).
if(exist([output_path_log,'counterfactual.log']))
    delete([output_path_log,'counterfactual.log']);
end
% Write new log-file.
diary([output_path_log,'counterfactual.log']);
fprintf('This log contains background information on the counterfactuals based on model %d:\n',model_type);
fprintf('Estimates from: %d \n 1 for efficient GMM estimates \n 0 for first stage 2SLS GMM estimates \n',efficient_GMM);
fprintf(strcat('Counterfactuals ran on: ',date,'\n'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Export observed and predicted firm-level market shares.
ms_obs = reshape(data(:,5),[J,T])';
% Aggregate to firm level in order to avoid confidentiality issues.
ms_firm_obs = [ms_obs(:,1) + ms_obs(:,2), ms_obs(:,3) + ms_obs(:,4), ...
                 ms_obs(:,5), ms_obs(:,6) + ms_obs(:,7),  ...
                 ms_obs(:,8) + ms_obs(:,9), ms_obs(:,10), ms_obs(:,11)];
% Export predicted firm-level market shares (to highlight quality of model fit).
dlmwrite(project_paths('OUT_ANALYSIS','ms_firm_obs.csv'),ms_firm_obs);
% Compute average market shares over time.
ms_avg_firm_obs = mean(ms_firm_obs(2:end,:),1);

% Aggregate to firm level in order to avoid confidentiality issues.
ms_pred_obs = reshape(pred_shares1d,[J,T])';
ms_firm_pred = [ms_pred_obs(:,1) + ms_pred_obs(:,2), ms_pred_obs(:,3) + ms_pred_obs(:,4), ...
                 ms_pred_obs(:,5), ms_pred_obs(:,6) + ms_pred_obs(:,7),  ...
                 ms_pred_obs(:,8) + ms_pred_obs(:,9), ms_pred_obs(:,10), ms_pred_obs(:,11)];
% Export predicted firm-level market shares (to highlight quality of model fit).
dlmwrite(project_paths('OUT_ANALYSIS','ms_firm_pred.csv'),ms_firm_pred);
% Compute average market shares over time.
ms_avg_firm_pred = mean(ms_firm_pred(2:end,:),1);
% Compute difference between predicted and observed market shares.
% Just another safety check: This should be zero!
ms_diff_predobs = ms_pred_obs - ms_obs;
fprintf('Average difference between predicted and observed market share: %6.6f.\n',mean(mean(ms_diff_predobs)));
%% Compute aggregate consumer surplus in status quo observed in the data.
% This is measured in unit of price coefficient (here: 100-EUR)
% This is not the correct way of aggregating because it ignores the
% joint distribution of lagged shares and PCW usage decision. Correct way
% of computing aggregate welfare is below.
CS_data_agg_old = pred_PCWAgg1d .* cSurplusPCW_obs + (1.0-pred_PCWAgg1d) .* cSurplusNoPCW_obs; 
% Compute consumer surplus in data separately for different age groups.
CS_i_young_data = cSurplusI_obs(:,preference_shocks(1,:)'==0);
CS_i_old_data = cSurplusI_obs(:,preference_shocks(1,:)'==1);
% Recall that consumer surplsu statistics from model file are already
% scaled by (negative of) price coefficient, i.e., this is in 100-EUR.
% Consumer surplus of the young in observed data (in 100 EUR).
CS_i_young_data = mean( CS_i_young_data,2);
% Consumer surplus of the old in observed data (in 100 EUR).
CS_i_old_data = mean( CS_i_old_data,2);
% Sanity check: average welfare over all consumers and see whether we get
% the same as from model file.
CS_i_all_data = mean(cSurplusI_obs,2);
% This is the correct way of computing aggregate welfare observed in data.
CS_data_agg = CS_i_all_data;
fprintf('Compare the two aggregate welfare statistics as a sanity check:\n %6.4f \n %6.4f\n',mean(CS_data_agg_old),mean(CS_data_agg));
fprintf('Compare consumer surplus for all, young, and old consumers: \n %6.4f \n %6.4f \n %6.4f \n',mean(CS_i_all_data), mean(CS_i_young_data),mean(CS_i_old_data));
% Combine welfare statistics for both types in one matrix.
CS_data_agg = [CS_i_all_data, CS_i_young_data, CS_i_old_data];
CS_data_yt = [...
    mean(CS_data_agg(1:12,:)); ...
    mean(CS_data_agg(13:24,:)); ...
    mean(CS_data_agg(25:36,:)); ...
    mean(CS_data_agg(37:48,:)); ...
    mean(CS_data_agg(48:end,:))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% START ACTUAL COUNTERFACTUALS HERE.
run_counterfactuals = 1;
run_estimation = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Counterfactual 1: Effects of switching subsidy
cf_subsidy_switch = 1;
% Set subsidy (in 100 EUR): We can experiment with other numbers here.
% Baseline model estimates are 23 EUR, so set it to 95% of that: 21 EUR
% Extended model estimates are 16 and 21 EUR, so set it to 95% of average: 20 EUR
if model_type==1
    subsidy_switch = 0.21;
elseif model_type==4
    subsidy_switch = 0.20;
elseif model_type==5
    subsidy_switch=0;
elseif model_type==6
    subsidy_switch=0.21;
end
    
% Resolve model for new market share distribution when switching subsidy is granted.
fprintf('Computing counterfactual when switching costs are reduced:\n');
[~,~,~,~,~,~,~, ms_subsidy_long, cr_subsidy,CS_subsidy_switch_PCW,CS_subsidy_switch_noPCW,pcw_use_subsidy_switch, ...
    ~,~,~,~,~,cSurplusI_subsidy] = gmm_function(theta_opt_cf);
% Reshape counterfactual shares into wide format.
ms_subsidy = reshape(ms_subsidy_long,[J,T])';
% Write market share distribution with switching subsidy to file.
% Aggregate to firm level in order to avoid confidentiality issues.
ms_firm_subsidy = [ms_subsidy(:,1) + ms_subsidy(:,2), ms_subsidy(:,3) + ms_subsidy(:,4), ...
                 ms_subsidy(:,5), ms_subsidy(:,6) + ms_subsidy(:,7),  ...
                 ms_subsidy(:,8) + ms_subsidy(:,9), ms_subsidy(:,10), ms_subsidy(:,11)];
% Export predicted firm-level market shares (to highlight quality of model fit).
dlmwrite(project_paths('OUT_ANALYSIS',['ms_firm_subsidy_' , file_suffix,'.csv']),ms_firm_subsidy);
% Compute average market shares over time.
ms_avg_firm_subsidy = mean(ms_firm_subsidy,1);
% Compute difference between observed and counterfactual market shares.
% One-dimensional long format.
ms_diff_switch_subsidy = ms_subsidy_long - pred_shares1d;
% Wide matrix format (firms in columns).
ms_diff_switch_subsidy_agg = ms_subsidy - ms_pred_obs;

% Some sanity checks to debug counterfactual code.
ms_subsidy_placebo_diff = ms_subsidy_long - pred_shares1d;
hist(ms_subsidy_placebo_diff)

%% Compute costs of switching subsidy scenario.
switch_subsidy_costs = cr_subsidy .* subsidy_switch;
% Counterfactual consumer surplus is measured in unit of price coefficient (here: 100 EUR).
% Compute consumer surplus in data separately for different age groups.
CS_i_young_subsidy_switch = cSurplusI_subsidy(:,preference_shocks(1,:)'==0);
CS_i_old_subsidy_switch = cSurplusI_subsidy(:,preference_shocks(1,:)'==1);
CS_subsidy_switch_all = mean(cSurplusI_subsidy,2);
CS_subsidy_switch_young = mean(CS_i_young_subsidy_switch,2);
CS_subsidy_switch_old = mean(CS_i_old_subsidy_switch,2);
% Combine in one matrix.
CS_subsidy_switch_agg = [CS_subsidy_switch_all, CS_subsidy_switch_young, CS_subsidy_switch_old];


% Net gain in CS, i.e. substract cost of switching subsidy.
% Compensating variation is transformed to EUR/month, i.e., multiply model output by 100.
CS_CV_subsidy_switch = 100 * (CS_subsidy_switch_agg - CS_data_agg); 
fprintf('Average increase in consumer surplus per month and consumer (in EUR) in switching subsidy counterfactual \n Averaged over all consumers: \n %6.2f \n Averaged over young consumers: \n %6.2f \n Averaged over old consumers: \n %6.2f\n', ...
    mean(CS_CV_subsidy_switch(:,1)), mean(CS_CV_subsidy_switch(:,2)), mean(CS_CV_subsidy_switch(:,3)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Counterfactual 2: Effect of eliminating limited information.
cf_subsidy_switch = 0;
cf_perfect_info = 1;

% Factor by which to reduce PCW search cost.
kappa_cut = 0.95;
% Resolve for new market share distribution when all consumers are (more)
% informed via the PCW.
fprintf('Computing counterfactual when search costs are reduced:\n');
[~,~,~,~,~,~,~,ms_perfect_info_long, cr_perfect_info,CS_perfect_info_PCW,CS_perfect_info_noPCW,pcw_use_perfect_info, ...
    ~,~,~,~,~,cSurplusI_perfect_info] = gmm_function(theta_opt_cf);
% Reshape counterfactual market shares into wide format.
ms_perfect_info = reshape(ms_perfect_info_long,[J,T])';

% Write market share distribution with perfect information to file.
% Aggregate to firm level in order to avoid confidentiality issues.
ms_firm_perfect_info = [ms_perfect_info(:,1) + ms_perfect_info(:,2), ms_perfect_info(:,3) + ms_perfect_info(:,4), ...
                 ms_perfect_info(:,5), ms_perfect_info(:,6) + ms_perfect_info(:,7),  ...
                 ms_perfect_info(:,8) + ms_perfect_info(:,9), ms_perfect_info(:,10), ms_perfect_info(:,11)];
% Export predicted firm-level market shares (to highlight quality of model fit).
dlmwrite(project_paths('OUT_ANALYSIS',['ms_firm_perfect_info_',file_suffix,'.csv']),ms_firm_perfect_info);


% Compute average market shares over time.
ms_avg_firm_perfect_info = mean(ms_firm_perfect_info,1);
% Compute difference between observed and counterfactual market shares.
% One-dimensional long format.
ms_diff_perfect_info = ms_perfect_info_long - pred_shares1d;
% Wide matrix format (firms in columns).
ms_diff_perfect_info_agg = ms_perfect_info - ms_pred_obs;

% Compute compensating variation, i.e. EUR-amount for which the average consumers would
% be indifferent between actual market and counterfactual market structure.
% Multiply by 100 because price coefficient is measured for 100 EUR.
%% Compute costs for subsidizing search costs.
perfect_info_costs = pcw_use_perfect_info .* kappa_cut .* kappa_out_EUR(:,1) ./ 100 ;

% Compute consumer surplus in data separately for different age groups.
CS_i_young_perfect_info = cSurplusI_perfect_info(:,preference_shocks(1,:)'==0);
CS_i_old_perfect_info = cSurplusI_perfect_info(:,preference_shocks(1,:)'==1);
CS_perfect_info_all = mean(cSurplusI_perfect_info,2);
CS_perfect_info_young = mean(CS_i_young_perfect_info,2);
CS_perfect_info_old = mean(CS_i_old_perfect_info,2);
% Combine in one matrix.
CS_perfect_info_agg = [CS_perfect_info_all, CS_perfect_info_young, CS_perfect_info_old];
% Compute compensatig variation.
CS_CV_perfect_info = 100 * (CS_perfect_info_agg - CS_data_agg); 
fprintf('Average increase in consumer surplus per month and consumer (in EUR) in perfect info counterfactual \n Averaged over all consumers: \n %6.2f \n Averaged over young consumers: \n %6.2f \n Averaged over old consumers: \n %6.2f\n', ...
    mean(CS_CV_perfect_info(:,1)), mean(CS_CV_perfect_info(:,2)), mean(CS_CV_perfect_info(:,3)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Counterfactual 3: Evaluate effect of regulator campaign.
% Set kappa parameters as if there was no regulator campaign that reduced kappa search costs. 
cf_subsidy_switch = 0;
cf_perfect_info = 0;
cf_reg_campaign = 1;
% Resolve for new market share distribution when all consumers are
% informed.
fprintf('Computing counterfactual when regulator campaign is eliminated:\n');
[~,~,~,~,~,~,~,ms_reg_campaign_long, cr_reg_campaign,CS_reg_campaign_PCW,CS_reg_campaign_noPCW,pcw_use_reg_campaign, ...
    ~,~,~,~,~, cSurplusI_campaign] = gmm_function(theta_opt_cf);
% Aggregate counterfactual market shares.
ms_reg_campaign = reshape(ms_reg_campaign_long,[J,T])';

% Write market share distribution with perfect information to file.
% Aggregate to firm level in order to avoid confidentiality issues.
ms_firm_reg_campaign = [ms_reg_campaign(:,1) + ms_reg_campaign(:,2), ms_reg_campaign(:,3) + ms_reg_campaign(:,4), ...
                 ms_reg_campaign(:,5), ms_reg_campaign(:,6) + ms_reg_campaign(:,7),  ...
                 ms_reg_campaign(:,8) + ms_reg_campaign(:,9), ms_reg_campaign(:,10), ms_reg_campaign(:,11)];
% Export predicted firm-level market shares (to highlight quality of model fit).
dlmwrite(project_paths('OUT_ANALYSIS','ms_firm_reg_campaign.csv'),ms_firm_reg_campaign);

% Compute average market shares over time.
ms_avg_firm_reg_campaign = mean(ms_firm_reg_campaign,1);
% Compute difference between observed and counterfactual market shares.
ms_diff_reg_campaign = ms_reg_campaign_long - pred_shares1d;
ms_diff_reg_campaign_agg = ms_reg_campaign - ms_pred_obs;
% Look at averages only during active campaign.
ms_avg_firm_reg_campaign_active = mean(ms_firm_reg_campaign(9:14,:));
ms_avg_firm_obs_campaign_active = mean(ms_firm_obs(9:14,:));
cf_reg_campaign = 0;

% Compute compensating variation, i.e. EUR-amount for which the average consumers would
% be indifferent between actual market and counterfactual market structure.
% Multiply by 100 because price coefficient is measured for 100 EUR.
%% Welfare gain from conducting regulator campaign.
% Compute consumer surplus in data separately for different age groups.
CS_i_young_campaign = cSurplusI_campaign(:,preference_shocks(1,:)'==0);
CS_i_old_campaign = cSurplusI_campaign(:,preference_shocks(1,:)'==1);
CS_campaign_all = mean(cSurplusI_campaign,2);
CS_campaign_young = mean(CS_i_young_campaign,2);
CS_campaign_old = mean(CS_i_old_campaign,2);
% Combine in one matrix.
CS_campaign_agg = [CS_campaign_all, CS_campaign_young, CS_campaign_old];
% Compute compensatig variation.
CS_CV_campaign = 100 * (CS_data_agg - CS_campaign_agg); 
fprintf('Average increase in consumer surplus per month and consumer (in EUR) in regulator campaign counterfactual \n Averaged over all consumers: \n %6.2f \n Averaged over young consumers: \n %6.2f \n Averaged over old consumers: \n %6.2f\n', ...
    mean(CS_CV_campaign(:,1)), mean(CS_CV_campaign(:,2)), mean(CS_CV_campaign(:,3)));
% Scale consumer surplus gain in EUR.
% CS_CV_reg_campaign = 100 * (CS_data_agg - CS_reg_campaign_agg);
% Print welfare gains form regulator campaign by year:
% Print welfare gains per average consumer
% Pre-campaign period: months 1 - 8 (should be zero)
% Active campaign period: months 9 - 14
% Six months right after campaign: months 15-20
% First year post-campaign: months 21-32
% Second year post-campaign: months 33-44
% Last months of sample: months 45-53
fprintf('Total welfare gain per consumer before campaign (in EUR):\n%6.2f, %6.2f, %6.2f \n',sum(CS_CV_campaign(1:7,:),1));
fprintf('Total welfare gain per consumer during active campaign (in EUR):\n%6.2f, %6.2f, %6.2f \n',sum(CS_CV_campaign(8:13,:),1));
fprintf('Total welfare gain per consumer six months after campaign (in EUR):\n%6.2f, %6.2f, %6.2f \n',sum(CS_CV_campaign(14:20,:),1));
fprintf('Total welfare gain per consumer first full year after campaign (in EUR):\n%6.2f, %6.2f, %6.2f \n',sum(CS_CV_campaign(21:32,:),1));
fprintf('Total welfare gain per consumer second full year after campaign (in EUR):\n%6.2f, %6.2f, %6.2f \n',sum(CS_CV_campaign(33:44,:),1));
fprintf('Total welfare gain per consumer during last 9 months of sample (in EUR):\n%6.2f, %6.2f, %6.2f \n',sum(CS_CV_campaign(45:end,:),1));
fprintf('Total welfare gain per consumer during whole sample period (in EUR):\n%6.2f, %6.2f, %6.2f \n',sum(CS_CV_campaign));
fprintf('Total welfare gain for all of Flanders before campaign (in mio EUR):\n%6.2f, %6.2f, %6.2f \n',2.7 * sum(CS_CV_campaign(1:7,:),1));
fprintf('Total welfare gain for all of Flanders during active campaign (in mio EUR):\n%6.2f, %6.2f, %6.2f \n',2.7 * sum(CS_CV_campaign(8:13,:),1));
fprintf('Total welfare gain for all of Flanders six months after campaign (in mio EUR):\n%6.2f, %6.2f, %6.2f \n',2.7 * sum(CS_CV_campaign(14:20,:),1));
fprintf('Total welfare gain for all of Flanders first full year after campaign (in mio EUR):\n%6.2f, %6.2f, %6.2f \n',2.7 * sum(CS_CV_campaign(21:32,:),1));
fprintf('Total welfare gain for all of Flanders second full year after campaign (in mio EUR):\n%6.2f, %6.2f, %6.2f \n',2.7 * sum(CS_CV_campaign(33:44,:),1));
fprintf('Total welfare gain for all of Flanders during last 9 months of sample (in mio EUR):\n%6.2f, %6.2f, %6.2f \n',2.7 * sum(CS_CV_campaign(45:end,:),1));
fprintf('Total welfare gain for all of Flanders during whole sample period (in mio EUR):\n%6.2f, %6.2f, %6.2f \n',2.7 * sum(CS_CV_campaign,1));

%% Overview of different market share distributions in different counterfactuals.
% 1. observed shares
% 2. predicted shares
% 3. shares with switching subsidy
% 4. shares with perfect information subsidy
% 5. shares with regulator campaign eliminated (averaged over whoel sample).
% 6. shares without regulator campaign (averaged only over months of
% campaign
% 7. observed average shares during regulator campaign.

ms_avg_firm_total = [ ...
    ms_avg_firm_obs; ...
    ms_avg_firm_pred; ...
    ms_avg_firm_subsidy; ...
    ms_avg_firm_perfect_info; ...
    ms_avg_firm_reg_campaign; ...
    ms_avg_firm_reg_campaign_active; ...
    ms_avg_firm_obs_campaign_active];
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Counterfactuals: What if no market liberalization? 
% Call model file to get consumer-specific utility mu and
% individual-specific price coefficients.
delta_opt = log(expmvalold);
% (A) First case: compute consumer surplus when only incumbent contract is available.
% Take existing deltas for incumbent.
delta_wide = reshape(delta_opt,J,T)';
delta_inc = repmat(delta_wide(:,1),1,NS_cons);
delta_inc_green = repmat(delta_wide(:,2),1,NS_cons);

% Extract incumbent mu.
mu_wide = reshape(mu_obs,[J,T,NS_cons]);
mu_inc = squeeze(mu_wide(1,:,:));
mu_inc_green = squeeze(mu_wide(2,:,:));
% Compute total utility form only having incumbent contract available.
v_inc_cf_1 = delta_inc + mu_inc;
v_inc_cf_1_green = delta_inc_green + mu_inc_green;
v_inc_cf_exp = exp(cat(3,v_inc_cf_1, v_inc_cf_1_green));

% Extract valuations for young and old types separately.
v_young_inc_only_1 = v_inc_cf_1(:,preference_shocks(1,:)'==0);
v_old_inc_only_1 = v_inc_cf_1(:,preference_shocks(1,:)'==1);
v_young_inc_only_1_green = v_inc_cf_1_green(:,preference_shocks(1,:)'==0);
v_old_inc_only_1_green = v_inc_cf_1_green(:,preference_shocks(1,:)'==1);
price_coeff_young = price_coeff_obs(preference_shocks(1,:)'==0);  
price_coeff_old = price_coeff_obs(preference_shocks(1,:)'==1);  

%% If switching costs are contract specific, we have to solve for actual market shares (because past contract affects current surplus)
if contract_sc==1 % case where switching costs work also within firms/across contracts.
    % When there is just one contract, no switching is possible so
    % switching costs are not relevant at all.
    % Compute consumer surplus in data separately for different age groups.
    % Since here there is only one product in the choice set log(exp())
    % just cancels out.
    CS_inc_only_1_all = - mean(v_inc_cf_1 ./ repmat(price_coeff_obs',T,1),2);
    CS_inc_only_1_young = - mean(v_young_inc_only_1 ./ repmat(price_coeff_young',T,1),2);
    CS_inc_only_1_old = - mean(v_old_inc_only_1 ./ repmat(price_coeff_old',T,1),2);
    % Combine in one matrix.
    CS_inc_only_1_agg = [CS_inc_only_1_all, CS_inc_only_1_young, CS_inc_only_1_old];

    % For robustness check where incumbent offers conventional contract 11 times.
    % Keep in mind that there are assumed to be switching costs between the different pseudo-contracts.
    % Define payoffs from logit-equivalent portfolio.
    v_inc_cf_exp_1_logitJ = repmat(v_inc_cf_exp(:,:,1),[1,1,J]);
    % Note that the CS matrix already contains three separate columns for averages over all, young, old consumers, respectively.
    [s_i_inc_only_1_logitJ, CS_inc_only_1_logitJ_agg] = ...
        compute_cf_shares(theta_opt_cf,preference_shocks,v_inc_cf_exp_1_logitJ,price_coeff_obs);
    
    % Compute welfare separately for each year and consumer type (young and
    % old).
    CS_inc_only_1_yt = [...
    mean(CS_inc_only_1_agg(1:12,:)); ...
    mean(CS_inc_only_1_agg(13:24,:)); ...
    mean(CS_inc_only_1_agg(25:36,:)); ...
    mean(CS_inc_only_1_agg(37:48,:)); ...
    mean(CS_inc_only_1_agg(48:end,:))];

    CS_inc_only_1_logitJ_yt = [...
    mean(CS_inc_only_1_logitJ_agg(1:12,:)); ...
    mean(CS_inc_only_1_logitJ_agg(13:24,:)); ...
    mean(CS_inc_only_1_logitJ_agg(25:36,:)); ...
    mean(CS_inc_only_1_logitJ_agg(37:48,:)); ...
    mean(CS_inc_only_1_logitJ_agg(48:end,:))];
    
    % When incumbent offers both conventional and green contract, switching
    % costs matter, model needs to be resolved explicitly.
    [s_i_inc_only_2, CS_inc_only_2_agg] = ...
        compute_cf_shares(theta_opt_cf,preference_shocks,v_inc_cf_exp,price_coeff_obs);
    % What happens when incumbent offers contract portfolio 11 times.
    % Define payoffs from logit-equivalent portfolio.
    v_inc_cf_exp_logitJ = cat(3,repmat(v_inc_cf_exp(:,:,1),[1,1,5]),repmat(v_inc_cf_exp(:,:,2),[1,1,6]));
    [s_i_inc_only_2_logitJ, CS_inc_only_2_logitJ_agg] = ...
        compute_cf_shares(theta_opt_cf,preference_shocks,v_inc_cf_exp_logitJ,price_coeff_obs);
    
    % What if I only have incumbent contracts available and I can switch
    % freely across the two contracts?
    % Counterfactual variation 2: Incumbent offers both conventional and green
    % contract and consumers can switch back and forth without any frictions.
    CS_inc_only_2_nosc_all = - mean(  log(exp(v_inc_cf_1) + exp(v_inc_cf_1_green)) ./ repmat(price_coeff_obs',T,1),2);
    CS_inc_only_2_nosc_young = - mean(  log(exp(v_young_inc_only_1) + exp(v_young_inc_only_1_green)) ./ repmat(price_coeff_young',T,1),2);
    CS_inc_only_2_nosc_old = - mean(  log(exp(v_old_inc_only_1) + exp(v_old_inc_only_1_green)) ./ repmat(price_coeff_old',T,1),2);
    % Combine in one matrix.
    CS_inc_only_2_nosc_agg = [CS_inc_only_2_nosc_all, CS_inc_only_2_nosc_young, CS_inc_only_2_nosc_old];
    CS_inc_only_2_nosc_yt = [...
    mean(CS_inc_only_2_nosc_agg(1:12,:)); ...
    mean(CS_inc_only_2_nosc_agg(13:24,:)); ...
    mean(CS_inc_only_2_nosc_agg(25:36,:)); ...
    mean(CS_inc_only_2_nosc_agg(37:48,:)); ...
    mean(CS_inc_only_2_nosc_agg(48:end,:))];
    CS_CV_inc_only_2_nosc = 100 * (CS_inc_only_2_nosc_agg - CS_data_agg);
    fprintf('Average change in consumer surplus per month and consumer (in EUR) \nwhen only the incumbent contracts (conventional and green) as in the data is offered\nwith no switching frictions: \n %6.2f,%6.2f,%6.2f \n',mean(CS_CV_inc_only_2_nosc));

    
    % Finally, do logit correction for a setting where also switching costs
    % are eliminated. This could be compared to best-case scenario in which
    % all contracts are offered without any frictions.
    
    CS_inc_only_2_logitJ_nosc_all = - mean(  log(5 * exp(v_inc_cf_1) + 6 * exp(v_inc_cf_1_green)) ./ repmat(price_coeff_obs',T,1),2);
    CS_inc_only_2_logitJ_nosc_young = - mean(  log(5*exp(v_young_inc_only_1) + 6* exp(v_young_inc_only_1_green)) ./ repmat(price_coeff_young',T,1),2);
    CS_inc_only_2_logitJ_nosc_old = - mean(  log(5*exp(v_old_inc_only_1) + 6*exp(v_old_inc_only_1_green)) ./ repmat(price_coeff_old',T,1),2);
    CS_inc_only_2_logitJ_nosc_agg = [ CS_inc_only_2_logitJ_nosc_all  CS_inc_only_2_logitJ_nosc_young CS_inc_only_2_logitJ_nosc_old];
    CS_CV_inc_only_2_logitJ_nosc = 100 * (CS_inc_only_2_logitJ_nosc_agg - CS_data_agg);
    fprintf('Average change in consumer surplus per month and consumer (in EUR) \nwhen only the incumbent contracts (conventional and green) as in the data is offered\nwith no switching frictions AND LOGIT CORRECTION APPLIED:\n \n %6.2f,%6.2f,%6.2f \n',mean(CS_CV_inc_only_2_logitJ_nosc));

    % Average separately for consumer types and years.
    CS_inc_only_2_yt = [...
    mean(CS_inc_only_2_agg(1:12,:)); ...
    mean(CS_inc_only_2_agg(13:24,:)); ...
    mean(CS_inc_only_2_agg(25:36,:)); ...
    mean(CS_inc_only_2_agg(37:48,:)); ...
    mean(CS_inc_only_2_agg(48:end,:))];
    % Controlling for logit shock equivalence.    
    CS_inc_only_2_logitJ_yt = [...
    mean(CS_inc_only_2_logitJ_agg(1:12,:)); ...
    mean(CS_inc_only_2_logitJ_agg(13:24,:)); ...
    mean(CS_inc_only_2_logitJ_agg(25:36,:)); ...
    mean(CS_inc_only_2_logitJ_agg(37:48,:)); ...
    mean(CS_inc_only_2_logitJ_agg(48:end,:))]; 
elseif contract_sc==0 % case where switching cost only work across firms.
    % Compute consumer surplus in data separately for different age groups.
    CS_inc_only_1_all = - mean(v_inc_cf_1 ./ repmat(price_coeff_obs',T,1),2);
    CS_inc_only_1_young = - mean(v_young_inc_only_1 ./ repmat(price_coeff_young',T,1),2);
    CS_inc_only_1_old = - mean(v_old_inc_only_1 ./ repmat(price_coeff_old',T,1),2);
    % Combine in one matrix.
    CS_inc_only_1_agg = [CS_inc_only_1_all, CS_inc_only_1_young, CS_inc_only_1_old];
    % Compute welfare separately for each year and consumer type (young and
    % old).
    CS_inc_only_1_yt = [...
    mean(CS_inc_only_1_agg(1:12,:)); ...
    mean(CS_inc_only_1_agg(13:24,:)); ...
    mean(CS_inc_only_1_agg(25:36,:)); ...
    mean(CS_inc_only_1_agg(37:48,:)); ...
    mean(CS_inc_only_1_agg(48:end,:))];

    % Counterfactual variation 2: Incumbent offers both conventional and green
    % contract and consumers can switch back and forth without any frictions.
    CS_inc_only_2_all = - mean(  log(exp(v_inc_cf_1) + exp(v_inc_cf_1_green)) ./ repmat(price_coeff_obs',T,1),2);
    CS_inc_only_2_young = - mean(  log(exp(v_young_inc_only_1) + exp(v_young_inc_only_1_green)) ./ repmat(price_coeff_young',T,1),2);
    CS_inc_only_2_old = - mean(  log(exp(v_old_inc_only_1) + exp(v_old_inc_only_1_green)) ./ repmat(price_coeff_old',T,1),2);
    % Combine in one matrix.
    CS_inc_only_2_agg = [CS_inc_only_2_all, CS_inc_only_2_young, CS_inc_only_2_old];
    
    CS_inc_only_2_yt = [...
    mean(CS_inc_only_2_agg(1:12,:)); ...
    mean(CS_inc_only_2_agg(13:24,:)); ...
    mean(CS_inc_only_2_agg(25:36,:)); ...
    mean(CS_inc_only_2_agg(37:48,:)); ...
    mean(CS_inc_only_2_agg(48:end,:))];
    
    
    % New counterfactual: bound welfare effects of logit shocks.
    % idea: incumbent offers identical contract 11 times so that consumers get
    % 11 different logit shocks. Here we assume that incumbent offers identical
    % contracts to observed market structure (i.e., 5 conventional and 6 green
    % contracts).
    % For robustness check where incumbent offers contract 11 times.
    CS_inc_only_1_logitJ = - mean( ...
        log( 5.* exp(v_inc_cf_1) + 6.* exp(v_inc_cf_1_green)) ...
        ./ repmat(price_coeff_obs',T,1),2);
    CS_inc_only_1_logitJ_young = - mean( ...
        log( 5.* exp(v_young_inc_only_1) + 6 .* exp(v_young_inc_only_1_green)) ...
        ./ repmat(price_coeff_young',T,1),2);
    CS_inc_only_1_logitJ_old = - mean( ...
        log( 5.* exp(v_old_inc_only_1) + 6.* exp(v_old_inc_only_1_green)) ...
        ./ repmat(price_coeff_old',T,1),2);
    % Combine in one matrix.
    CS_inc_only_1_logitJ_agg = [CS_inc_only_1_logitJ_all, CS_inc_only_1_logitJ_young, CS_inc_only_1_logitJ_old];
    
    CS_inc_only_1_logitJ_yt = [...
    mean(CS_inc_only_1_logitJ_agg(1:12,:)); ...
    mean(CS_inc_only_1_logitJ_agg(13:24,:)); ...
    mean(CS_inc_only_1_logitJ_agg(25:36,:)); ...
    mean(CS_inc_only_1_logitJ_agg(37:48,:)); ...
    mean(CS_inc_only_1_logitJ_agg(48:end,:))];   
end

% Compute compensating variation and print results.
CS_CV_inc_only_1 = 100 * (CS_inc_only_1_agg - CS_data_agg);
fprintf('Average change in consumer surplus per month and consumer (in EUR) \nwhen only the conventional incumbent contract as in the data is offered: \n %6.2f, %6.2f, %6.2f  \n',mean(CS_CV_inc_only_1));
CS_CV_inc_only_2 = 100 * (CS_inc_only_2_agg - CS_data_agg);
fprintf('Average change in consumer surplus per month and consumer (in EUR) \nwhen only the incumbent contracts (conventional and green) as in the data is offered: \n %6.2f,%6.2f,%6.2f \n',mean(CS_CV_inc_only_2));
CS_CV_inc_only_1_logitJ = 100 * (CS_inc_only_1_logitJ_agg - CS_data_agg);
fprintf('Average change in consumer surplus per month and consumer (in EUR) \nwhen incumbent offers identical contract J times: \n %6.2f,%6.2f,%6.2f \n',mean(CS_CV_inc_only_1_logitJ));
CS_CV_inc_only_2_logitJ = 100 * (CS_inc_only_2_logitJ_agg - CS_data_agg);
fprintf('Average change in consumer surplus per month and consumer (in EUR) \nwhen incumbent offers conventional and green contract J times: \n %6.2f,%6.2f,%6.2f \n',mean(CS_CV_inc_only_2_logitJ));

% Compute compensating variation for year and type bins.
CS_CV_inc_only_1_yt = 100*(CS_inc_only_1_yt - CS_data_yt);
CS_CV_inc_only_2_yt = 100*(CS_inc_only_2_yt - CS_data_yt);
CS_CV_inc_only_1_logitJ_yt = 100*(CS_inc_only_1_logitJ_yt - CS_data_yt);
CS_CV_inc_only_2_logitJ_yt = 100*(CS_inc_only_2_logitJ_yt - CS_data_yt);

% (B) Second case: compute consumer surplus when regulated utility resells wholesale + fixed markup.
% Compute new price charged by consumer. 
% Markup that firm charges to consumers (in percent above wholesale spot
% price).
% Load correctd wholesale price time series.
wh_data = readmatrix(project_paths('OUT_DATA','wholesale_corr_export.csv'));
% Select correct wholesale price measure based on ACER report and scale to
% monthly average consumption of household.
whp_acer = wh_data(2:end,3) .* 3500 ./ 1000 ./ 12 ./ 100;

% Compute welfare for grid of markups.
markup_grid = linspace(0,250,751)' ./100;
CS_CV_inc_only_1_WHP_grid = zeros(T,3,length(markup_grid));

% Some variables needed for consumer surplus calculations.
% Do we have to add vat here?
% Potentially: incorporate temporary VAT cut in 2014/15
%vat_rate = 0.21;
% Detailed VAT rates.
% From April 2014 (corrsponds to month 27) to August 2015 (corrsponds to month 44): temporary decrease from 21% to 6%.
% This will be reflected in our hypothetical retail prices.
% Not sure if this is what we should compare the situation to or whether we
% should assum VAT cut to be special to observed data that would not have
% materialized in counterfactual.
vat_rate = 0.21 * ones(T,1);
vat_rate(27:44,1) = 0.06;
vat_rate_wide = (1.0+repmat(vat_rate,1,J));
% Transform wholesale price in data matrix to be consistent with retail
% price data, i.e. convert mwh to kwh, scale to average consumption level
% of 3500kwh per year, and transform into monthly price in 100-EUR, finally
% add VAT.
whprice_aux = reshape(data(:,9),J,T)';
p_wh_old = (1.0 + vat_rate) .* (whprice_aux .* 3500 ./ 1000 ./ 12 ./ 100);

p_wh = repmat(whp_acer,1,J);

price_ret_wide = reshape(data(:,7),J,T)./100;
% Compute observed markups.
markup_ret = ((price_ret_wide' ./vat_rate_wide) - p_wh) ./ (price_ret_wide'./vat_rate_wide);
markup_wh = ((price_ret_wide' ./vat_rate_wide) - p_wh) ./ p_wh;
hist(markup_ret)
title('Distribution of markups relative to retail price')
hist(markup_wh)
title('Distribution of markups relative to wholesale price')
% Compute average markup charged by firms in the data.
mu_obs_avg = mean(markup_wh,1);

% For sanity check: Compare observed retail prices to wholesale spot price in every month.
% Need to add outside good price to price in data (which is already
% normalized for outside good).
price_data_check = reshape(X_reg(:,2),J,T)';
price_retMarkup_grid = zeros(T,J,length(markup_grid));

% Initialize container for type-specific consumer surplsu with WHP
% structure.
CS_inc_only_1_WHP = zeros(T,3);
% Compute actual change in consumer surplus for a given markup level.
for mu_idx=1:length(markup_grid)
    markup = markup_grid(mu_idx);
    % Add markup to correctly scaled wholesale price measure and add VAT.
    p_whmu = p_wh .* (1.0 + markup) .* vat_rate_wide;
    p_whmu_norm = p_whmu - repmat(p_og,1,J);
    % Compute actual markups charged by firms (retail price over wholesale
    % price).
    price_retMarkup_grid(:,:,mu_idx) = (price_data_check - p_whmu_norm) ./ (price_data_check + repmat(p_og,1,J));
    % Compute consumer surplus by adding price differential between observed and
    % hypothetical incumbent price.
    price_change = p_whmu_norm(:,1) - price_data_check(:,1);
    v_inc_cf_1_WHP_all = v_inc_cf_1 +  (repmat(price_change,1,NS_cons) .* repmat(price_coeff_obs',T,1));
    v_inc_cf_1_WHP_young = v_inc_cf_1_WHP_all(:,preference_shocks(1,:)'==0);
    v_inc_cf_1_WHP_old = v_inc_cf_1_WHP_all(:,preference_shocks(1,:)'==1);
    
    % Compute new consumer surplus for different consumer types.
    % Col 1: Averaged over all consumers
    % Col 2: Averaged only over young consumers.
    % Col 3: Averaged only over old consumers.
    CS_inc_only_1_WHP(:,1) = - mean((v_inc_cf_1_WHP_all) ./ repmat(price_coeff_obs',T,1),2);
    CS_inc_only_1_WHP(:,2) = - mean((v_inc_cf_1_WHP_young) ./ repmat(price_coeff_young',T,1),2);
    CS_inc_only_1_WHP(:,3) = - mean((v_inc_cf_1_WHP_old) ./ repmat(price_coeff_old',T,1),2);
    CS_CV_inc_only_1_WHP = 100 * (CS_inc_only_1_WHP - CS_data_agg);
    CS_CV_inc_only_1_WHP_grid(:,:,mu_idx) = CS_CV_inc_only_1_WHP;   
end

% Plot size of consumer surplus change for different markup levels and consumer types.
subplot(1,3,1)
CS_CV_inc_only_1_WHP_plot = squeeze(mean(CS_CV_inc_only_1_WHP_grid,1));
plot(markup_grid, CS_CV_inc_only_1_WHP_plot(1,:));
title('All consumers');
sgtitle('Change in consumer surplus for wholesale-plus contract (incumbent only)')
ylabel('Avg change in consumer surplus per month (in EUR)')
xlabel('Markup charged on wholesale spot price')
subplot(1,3,2)
plot(markup_grid, CS_CV_inc_only_1_WHP_plot(2,:));
title('Young consumers');
ylabel('Avg change in consumer surplus per month (in EUR)')
xlabel('Markup charged on wholesale spot price')
subplot(1,3,3)
plot(markup_grid, CS_CV_inc_only_1_WHP_plot(3,:));
title('Old consumers');
sgtitle('Change in consumer surplus for wholesale-plus contract (incumbent only)')
ylabel('Avg change in consumer surplus per month (in EUR)')
xlabel('Markup charged on wholesale spot price')
                 
% Find consumer surplus change closest to zero.
[CS_CV_inc_only_1_WHP_sol, CS_CV_inc_only_1_WHP_sol_idx] = min(abs(CS_CV_inc_only_1_WHP_plot),[],2);
CS_CV_inc_only_1_WHP_sol_mu = markup_grid(CS_CV_inc_only_1_WHP_sol_idx);
% Do actual markup comparison (of potential interest).
%price_MUChange = price_retMarkup_grid(:,:,CS_CV_inc_only_1_WHP_sol_idx);
fprintf('Consumers are indifferent between status quo and hypothetical monopolist\n with one contract when markup is %4.2f percent.\n',100*CS_CV_inc_only_1_WHP_sol_mu);
fprintf('In observed data, incumbent charges on average a markup (p-w)/w of %4.2f percent.\n',100*mu_obs_avg(1));

% Observed price time-series (in EUR) for incumbent contract.
price_cf_inc_1_WHP_only = [price_ret_wide(1,:)' * 100, 100 * (1.0+vat_rate) .*  p_wh(:,1) .* (1.0 + CS_CV_inc_only_1_WHP_sol_mu(1))]; 
% Write to table for formatting in Python.
price_labels = {'Observed','Indiff. Markup'};
cf_inc_only1_id_price = table(price_cf_inc_1_WHP_only(:,1),price_cf_inc_1_WHP_only(:,2), 'VariableNames', price_labels);
writetable(cf_inc_only1_id_price, project_paths('OUT_ANALYSIS',['cf_inc_only1_price_id_',file_suffix_full,'.csv']));

close all
subplot(1,1,1)
plot(time_grid, price_cf_inc_1_WHP_only(:,1), time_grid, price_cf_inc_1_WHP_only(:,2));
title('Change in consumer surplus for wholesale-plus contract (incumbent only)')
ylabel('Monthly retail price (in EUR)')
xlabel('Time')
legend('Observed','Indifference markup')
saveas(gcf,project_paths('OUT_FIGURES',['prices_cf_inc_1_WHP_',file_suffix,'.pdf']))


%% New counterfactual: compute indifference markup when incumbent offers two contracts.
% Compute new price charged by consumer. 
% Markup that firm charges to consumers (in percent above wholesale spot
% price).
% Compute welfare for grid of markups.
% Because here model has to be resolved, let's focus on fewer markup levels for exploration purposes.
markup_grid_2 = linspace(0,250,751)' ./100;
CS_CV_inc_2_WHP_grid = zeros(T,3,length(markup_grid_2));
% Detailed VAT rates that incorporate temporary VAT cut.
% From April 2014 (corrsponds to month 27) to August 2015 (corrsponds to month 44): temporary decrease from 21% to 6%.
% This will be reflected in our hypothetical retail prices.
vat_rate = 0.21 * ones(T,1);
vat_rate(27:44,1) = 0.06;
% Transform wholesale price in data matrix to be consistent with retail
% price data, i.e. convert mwh to kwh, scale to average consumption level
% of 3500kwh per year, and transform into monthly price in 100-EUR, finally
% add VAT.
whprice_aux = reshape(data(:,9),J,T)';
%p_wh = (1.0 + vat_rate) .* (whprice_aux .* 3500 ./ 1000 ./ 12 ./ 100);
% Transform wholesale price in data matrix to be consistent with retail
% price data, i.e. convert mwh to kwh, scale to average consumption level
% of 3500kwh per year, and transform into monthly price in 100-EUR, finally
% add VAT.
%p_wh = (1.0 + vat_rate) .* (whprice_aux .* 3500 ./ 1000 ./ 12 ./ 100);
% TO ADD: Green energy surcharge for last counterfactual.
% SW: For now I assume that there is a 3 cent surcharge for each kwh of
% green electricity. We need to supply a corect number instead of 0.03 here.
% Then compute: 3500 kwh yearly consumption, scale in 100-EUR and divide by
% 12 to get monthly charge.

green_surcharge = green_surcharge_factor / 100 * 3500 / 12;
green_idx = repmat([0,1,0,1,1,0,1,0,1,1,0],T,1);
p_wh = p_wh + green_surcharge .* green_idx;
% Normalize for outside good price.
p_og_long = kron(p_og,ones(J,1));
og_idx = repmat([zeros(J-1,1);1],T,1);
% Effective, i.e., normalized by outside good price, needs to be in long-format because model is resolved below.
p_wh_long = reshape(p_wh',J*T,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute actual change in consumer surplus for a given markup level.
for mu_idx=1:length(markup_grid_2)
    markup = markup_grid_2(mu_idx);
    % Add markup to correctly scaled wholesale price measure.
    p_whmu = p_wh .* (1.0 + markup) .* vat_rate_wide;
    p_whmu_norm = p_whmu - repmat(p_og,1,J);
    % Compute consumer surplus by adding price differential between observed and
    % hypothetical incumbent price.
    price_change = p_whmu_norm - price_data_check;
    price_util_adjustment = (repmat(price_change,1,1,NS_cons) .* permute(repmat(price_coeff_obs,1,T,J),[2,3,1]));
        
    v_inc_cf_loop_aux = repmat(delta_wide,[1,1,NS_cons]) + permute(mu_wide,[2,1,3]) + price_util_adjustment;
    v_inc_cf_2 = permute(exp(v_inc_cf_loop_aux(:,1:2,:)),[1,3,2]);
    
    % Now need to solve market share predictions, because consumers have
    % choice between conventional and green contract.
    [s_i_inc_only_2, CS_inc_only_2_WHP] = ...
        compute_cf_shares(theta_opt_cf,preference_shocks,v_inc_cf_2,price_coeff_obs);
    CS_CV_inc_2_WHP = 100 * (CS_inc_only_2_WHP - CS_data_agg);
    CS_CV_inc_2_WHP_grid(:,:,mu_idx) = CS_CV_inc_2_WHP;   
end

% Plot size of consumer surplus change for different markup levels and consumer types.
subplot(1,3,1)
CS_CV_inc_only_2_WHP_plot = squeeze(mean(CS_CV_inc_2_WHP_grid,1));
plot(markup_grid_2, CS_CV_inc_only_2_WHP_plot(1,:));
title('All consumers');
sgtitle('Change in consumer surplus for wholesale-plus contract (incumbent conv & green)')
ylabel('Avg change in consumer surplus per month (in EUR)')
xlabel('Markup charged on wholesale spot price')
subplot(1,3,2)
plot(markup_grid_2, CS_CV_inc_only_2_WHP_plot(2,:));
title('Young consumers');
ylabel('Avg change in consumer surplus per month (in EUR)')
xlabel('Markup charged on wholesale spot price')
subplot(1,3,3)
plot(markup_grid_2, CS_CV_inc_only_2_WHP_plot(3,:));
title('Old consumers');
ylabel('Avg change in consumer surplus per month (in EUR)')
xlabel('Markup charged on wholesale spot price')
% Find consumer surplus change closest to zero.
[CS_CV_inc_only_2_WHP_sol, CS_CV_inc_only_2_WHP_sol_idx] = min(abs(CS_CV_inc_only_2_WHP_plot),[],2);
CS_CV_inc_only_2_WHP_sol_mu = markup_grid(CS_CV_inc_only_2_WHP_sol_idx);

% Do actual markup comparison (of potential interest).
price_MUChange_2 = price_retMarkup_grid(:,:,CS_CV_inc_only_2_WHP_sol_idx);
fprintf('Consumers are indifferent between status quo and hypothetical monopolist\n with conv and green contract when markup is %4.2f percent.\n',100*CS_CV_inc_only_2_WHP_sol_mu);
% Observed price time-series (in EUR) for incumbent contract.
price_cf_inc_2_WHP_only = [price_ret_wide(1,:)' .* 100, 100 *p_wh(:,1) .* (1.0+vat_rate) .* (1.0 + CS_CV_inc_only_2_WHP_sol_mu(1))]; 
% Write to table for formatting in Python.
price_labels = {'Observed','Indiff. Markup'};
cf_inc_only_2_id_price = table(price_cf_inc_2_WHP_only(:,1),price_cf_inc_2_WHP_only(:,2), 'VariableNames', price_labels);
writetable(cf_inc_only_2_id_price, project_paths('OUT_ANALYSIS',['cf_inc_only_2_price_id_',file_suffix_full,'.csv']));
close all
subplot(1,1,1)
plot(time_grid, price_cf_inc_2_WHP_only(:,1), time_grid, price_cf_inc_2_WHP_only(:,2));
title('Change in consumer surplus for wholesale-plus contract (incumbent conv & green)')
ylabel('Monthly retail price (in EUR)')
xlabel('Time')
legend('Observed','Indifference markup')
saveas(gcf,project_paths('OUT_FIGURES',['prices_cf_inc_2_WHP_',file_suffix,'.pdf']))

%% Compute indifference markup when incumbent offers two contracts but no
% switching costs are paid.
%% Repeat exercise with various markups on wholesale price (wholesale-plus contracts)
% Compute new price charged by consumer. 
% Markup that firm charges to consumers (in percent above wholesale spot
% price).
% Compute welfare for grid of markups.
% Because here model has to be resolved, let's focus on fewer markup levels for exploration purposes.
markup_grid_2_nosc = linspace(000,250,751)' ./100;
CS_CV_inc_2_nosc_grid = zeros(T,3,length(markup_grid_2_nosc));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute actual change in consumer surplus for a given markup level.
CS_inc_2_nosc_WHP = zeros(T,3);
for mu_idx=1:length(markup_grid_2_nosc)
    markup = markup_grid_2_nosc(mu_idx);
    % Add markup to correctly scaled wholesale price measure.
    p_whmu = p_wh .* (1.0 + markup) .* vat_rate_wide;
    p_whmu_norm = p_whmu - repmat(p_og,1,J);
    % Compute consumer surplus by adding price differential between observed and
    % hypothetical incumbent price.
    price_change = p_whmu_norm - price_data_check;
    price_util_adjustment = (repmat(price_change,1,1,NS_cons) .* permute(repmat(price_coeff_obs,1,T,J),[2,3,1]));
    v_inc_cf_2_nosc_aux = repmat(delta_wide,[1,1,NS_cons]) + permute(mu_wide,[2,1,3]) + price_util_adjustment;
    % Keep only incumbent contracts.
    v_inc_cf_2_nosc = v_inc_cf_2_nosc_aux(:,1:2,:);
       
    v_inc_cf_2_nosc_WHP_all = squeeze(log(sum(exp(v_inc_cf_2_nosc),2)));
    v_inc_cf_2_nosc_WHP_young = v_inc_cf_2_nosc_WHP_all(:,preference_shocks(1,:)'==0);
    v_inc_cf_2_nosc_WHP_old = v_inc_cf_2_nosc_WHP_all(:,preference_shocks(1,:)'==1);

    % Compute new consumer surplus for different consumer types.
    % Col 1: Averaged over all consumers
    % Col 2: Averaged only over young consumers.
    % Col 3: Averaged only over old consumers.
    CS_inc_2_nosc_WHP(:,1) = - mean(v_inc_cf_2_nosc_WHP_all ./ repmat(price_coeff_obs',T,1),2);
    CS_inc_2_nosc_WHP(:,2) = - mean(v_inc_cf_2_nosc_WHP_young ./ repmat(price_coeff_young',T,1),2);
    CS_inc_2_nosc_WHP(:,3) = - mean(v_inc_cf_2_nosc_WHP_old ./ repmat(price_coeff_old',T,1),2);
    
    CS_CV_inc_2_nosc_WHP = 100 * (CS_inc_2_nosc_WHP - CS_data_agg);
    CS_CV_inc_2_nosc_grid(:,:,mu_idx) = CS_CV_inc_2_nosc_WHP;   
end
% Plot size of consumer surplus change for different markup levels and consumer types.
subplot(1,3,1)
sgtitle('Change in consumer surplus for wholesale-plus contract (incumbent conv & green, no switching costs)')
CS_CV_inc_2_nosc_plot = squeeze(mean(CS_CV_inc_2_nosc_grid,1));
plot(markup_grid_2_nosc, CS_CV_inc_2_nosc_plot(1,:));
title('All consumers');
ylabel('Avg change in consumer surplus per month (in EUR)')
xlabel('Markup charged on wholesale spot price')
subplot(1,3,2)
plot(markup_grid_2_nosc, CS_CV_inc_2_nosc_plot(2,:));
title('Young consumers');
ylabel('Avg change in consumer surplus per month (in EUR)')
xlabel('Markup charged on wholesale spot price')
subplot(1,3,3)
plot(markup_grid_2_nosc, CS_CV_inc_2_nosc_plot(3,:));
title('Old consumers');
ylabel('Avg change in consumer surplus per month (in EUR)')
xlabel('Markup charged on wholesale spot price')
% Find consumer surplus change closest to zero.
[CS_CV_inc_2_nosc_WHP_sol, CS_CV_inc_2_nosc_WHP_sol_idx] = min(abs(CS_CV_inc_2_nosc_plot),[],2);
CS_CV_inc_2_nosc_WHP_sol_mu = markup_grid(CS_CV_inc_2_nosc_WHP_sol_idx);



% Do actual markup comparison (of potential interest).
%price_MUChange = price_retMarkup_grid(:,:,CS_CV_inc_2_nosc_WHP_sol_idx);
fprintf('Consumers are indifferent between status quo and hypothetical monopolist\n with conv and green contracts (no switching costs) when markup is %4.2f - %4.2f - %4.2f percent.\n',100*CS_CV_inc_2_nosc_WHP_sol_mu);
% Observed price time-series (in EUR) for incumbent contract.
price_cf_inc_2_nosc_WHP = [price_ret_wide(1,:)' .* 100, 100 *p_wh(:,1) .* (1.0+vat_rate) .* (1.0 + CS_CV_inc_2_nosc_WHP_sol_mu(1))]; 
% Save raw price time series data for formattign in Python.
% Write to table for formatting in Python.
price_labels = {'Observed','Indiff. Markup'};
cf_inc_2_nosc_id_price = table(price_cf_inc_2_nosc_WHP(:,1),price_cf_inc_2_nosc_WHP(:,2), 'VariableNames', price_labels);
writetable(cf_inc_2_nosc_id_price, project_paths('OUT_ANALYSIS',['cf_inc_2_nosc_price_id_',file_suffix_full,'.csv']));
% Exploratory graph.
close all
subplot(1,1,1)
plot(time_grid, price_cf_inc_2_nosc_WHP(:,1), time_grid, price_cf_inc_2_nosc_WHP(:,2));
title('Change in consumer surplus for wholesale-plus contract (incumbent conv and green, no switching costs)')
ylabel('Monthly retail price (in EUR)')
xlabel('Time')
legend('Observed','Indifference markup')
saveas(gcf,project_paths('OUT_FIGURES',['prices_cf_inc_2_nosc_WHP_',file_suffix,'.pdf']))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute indifference markup when incumbent offers two contracts but no
% switching costs are paid and logit correction is applied.
%% Repeat exercise with various markups on wholesale price (wholesale-plus contracts)
% Compute new price charged by consumer. 
% Markup that firm charges to consumers (in percent above wholesale spot
% price).
% Compute welfare for grid of markups.
% Because here model has to be resolved, let's focus on fewer markup levels for exploration purposes.
markup_grid_2_lc_nosc = linspace(000,250,751)' ./100;
CS_CV_inc_2_lc_nosc_grid = zeros(T,3,length(markup_grid_2_lc_nosc));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute actual change in consumer surplus for a given markup level.
CS_inc_2_lc_nosc_WHP = zeros(T,3);
for mu_idx=1:length(markup_grid_2_lc_nosc)
    markup = markup_grid_2_lc_nosc(mu_idx);
    % Add markup to correctly scaled wholesale price measure.
    p_whmu = p_wh .* (1.0 + markup) .* vat_rate_wide;
    p_whmu_norm = p_whmu - repmat(p_og,1,J);
    % Compute consumer surplus by adding price differential between observed and
    % hypothetical incumbent price.
    price_change = p_whmu_norm - price_data_check;
    price_util_adjustment = (repmat(price_change,1,1,NS_cons) .* permute(repmat(price_coeff_obs,1,T,J),[2,3,1]));
    v_inc_cf_2_lc_nosc_aux = repmat(delta_wide,[1,1,NS_cons]) + permute(mu_wide,[2,1,3]) + price_util_adjustment;
    % Keep only incumbent contracts.
    v_inc_cf_2_lc_nosc = v_inc_cf_2_lc_nosc_aux(:,1:2,:);
    
  
    % v_inc_cf_2_lc_nosc_WHP_all = squeeze(log(sum(exp(v_inc_cf_2_nosc),2)));
    v_inc_cf_2_lc_nosc_WHP_all = squeeze(log( 5.* exp(v_inc_cf_2_lc_nosc(:,1,:)) + 6 .* exp(v_inc_cf_2_lc_nosc(:,2,:))));
    v_inc_cf_2_lc_nosc_WHP_young = v_inc_cf_2_lc_nosc_WHP_all(:,preference_shocks(1,:)'==0);
    v_inc_cf_2_lc_nosc_WHP_old = v_inc_cf_2_lc_nosc_WHP_all(:,preference_shocks(1,:)'==1);

    % Compute new consumer surplus for different consumer types.
    % Col 1: Averaged over all consumers
    % Col 2: Averaged only over young consumers.
    % Col 3: Averaged only over old consumers.
    CS_inc_2_lc_nosc_WHP(:,1) = - mean(v_inc_cf_2_lc_nosc_WHP_all ./ repmat(price_coeff_obs',T,1),2);
    CS_inc_2_lc_nosc_WHP(:,2) = - mean(v_inc_cf_2_lc_nosc_WHP_young ./ repmat(price_coeff_young',T,1),2);
    CS_inc_2_lc_nosc_WHP(:,3) = - mean(v_inc_cf_2_lc_nosc_WHP_old ./ repmat(price_coeff_old',T,1),2);
    
    CS_CV_inc_2_lc_nosc_WHP = 100 * (CS_inc_2_lc_nosc_WHP - CS_data_agg);
    CS_CV_inc_2_lc_nosc_grid(:,:,mu_idx) = CS_CV_inc_2_lc_nosc_WHP;   
end
% Plot size of consumer surplus change for different markup levels and consumer types.
subplot(1,3,1)
sgtitle('Change in consumer surplus for wholesale-plus contract (incumbent, logit corrected, no switching costs)')
CS_CV_inc_2_lc_nosc_plot = squeeze(mean(CS_CV_inc_2_lc_nosc_grid,1));
plot(markup_grid_2_lc_nosc, CS_CV_inc_2_lc_nosc_plot(1,:));
title('All consumers');
ylabel('Avg change in consumer surplus per month (in EUR)')
xlabel('Markup charged on wholesale spot price')
subplot(1,3,2)
plot(markup_grid_2_lc_nosc, CS_CV_inc_2_lc_nosc_plot(2,:));
title('Young consumers');
ylabel('Avg change in consumer surplus per month (in EUR)')
xlabel('Markup charged on wholesale spot price')
subplot(1,3,3)
plot(markup_grid_2_lc_nosc, CS_CV_inc_2_lc_nosc_plot(3,:));
title('Old consumers');
ylabel('Avg change in consumer surplus per month (in EUR)')
xlabel('Markup charged on wholesale spot price')
% Find consumer surplus change closest to zero.
[CS_CV_inc_2_lc_nosc_WHP_sol, CS_CV_inc_2_lc_nosc_WHP_sol_idx] = min(abs(CS_CV_inc_2_lc_nosc_plot),[],2);
CS_CV_inc_2_lc_nosc_WHP_sol_mu = markup_grid(CS_CV_inc_2_lc_nosc_WHP_sol_idx);

% Do actual markup comparison (of potential interest).
%price_MUChange = price_retMarkup_grid(:,:,CS_CV_inc_2_nosc_WHP_sol_idx);
fprintf('Consumers are indifferent between status quo and hypothetical monopolist\n with conv and green contracts (logit correction applied, no switching costs) when markup is %4.2f - %4.2f - %4.2f percent.\n',100*CS_CV_inc_2_lc_nosc_WHP_sol_mu);
% Observed price time-series (in EUR) for incumbent contract.
price_cf_inc_2_lc_nosc_WHP = [price_ret_wide(1,:)' .* 100, 100 *p_wh(:,1) .* (1.0+vat_rate) .* (1.0 + CS_CV_inc_2_lc_nosc_WHP_sol_mu(1))]; 
% Save raw price time series data for formattign in Python.
% Write to table for formatting in Python.
price_labels = {'Observed','Indiff. Markup'};
cf_inc_2_lc_nosc_id_price = table(price_cf_inc_2_lc_nosc_WHP(:,1),price_cf_inc_2_lc_nosc_WHP(:,2), 'VariableNames', price_labels);
writetable(cf_inc_2_lc_nosc_id_price, project_paths('OUT_ANALYSIS',['cf_inc_2_lc_nosc_price_id_',file_suffix_full,'.csv']));
% Exploratory graph.
close all
subplot(1,1,1)
plot(time_grid, price_cf_inc_2_lc_nosc_WHP(:,1), time_grid, price_cf_inc_2_lc_nosc_WHP(:,2));
title('Change in consumer surplus for wholesale-plus contract (incumbent conv and green, logit correction, no switching costs)')
ylabel('Monthly retail price (in EUR)')
xlabel('Time')
legend('Observed','Indifference markup')
saveas(gcf,project_paths('OUT_FIGURES',['prices_cf_inc_2_lc_nosc_WHP_',file_suffix,'.pdf']))

% End of indifference markup for logit corrected counterfactual in which
% incumbent offers identical copies of his own two contracts.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% New counterfactual: What if incumbent offers all contracts and gets rid of all market frictions?
% Incumbent offers all contracts at observed prices and eliminates all market frictions.
% This gives us an upper bound on what can be achieved in terms of welfare
% gain for consumers? This can be computed in closed form from estimates and doesn't require
% resolving the model.
% Compute observed utilities for each contract.
v_ij_all = log(squeeze(sum(reshape( ...
    exp(repmat(delta_opt,1,NS_cons) + mu_obs) ...
    ,[J,T,NS_cons]),1)));
v_ij_young = v_ij_all(:,preference_shocks(1,:)==0);
v_ij_old = v_ij_all(:,preference_shocks(1,:)==1);
% Compute consumer surplus for different consumer types (in unit of price
% coefficient, ie 100 EUR).
CS_inc_all_all = - mean(  v_ij_all  ./ repmat(price_coeff_obs',T,1),2);
CS_inc_all_young = - mean(  v_ij_young  ./ repmat(price_coeff_young',T,1),2);
CS_inc_all_old = - mean(  v_ij_old  ./ repmat(price_coeff_old',T,1),2);
CS_inc_all_agg= [CS_inc_all_all,CS_inc_all_young,CS_inc_all_old];

CS_CV_inc_all = 100 * (CS_inc_all_agg - CS_data_agg);
fprintf('Consumer surplus per month and consumer (in EUR) when incumbent offers all contracts without any  market frictions: \n All \t Young \t Old \n %4.2f \t %4.2f \t %4.2f \n',mean(CS_CV_inc_all));

% Compute market share predictions for case where incumbent ofers all
% contracts at observed prices.
share_inc_all = ... numerator
    reshape(exp(repmat(delta_opt,1,NS_cons) + mu_obs) ...
    ,[J,T,NS_cons]) ...
./ ... denominator
repmat(sum(reshape( ...
    exp(repmat(delta_opt,1,NS_cons) + mu_obs) ...
    ,[J,T,NS_cons]),1),[J,1,1]);

ms_avg_inc_all = squeeze(mean(mean(share_inc_all,3),2));
ms_avg_firm_inc_all =  ... % Aggregate to firm level in order to avoid confidentiality issues.
                [ms_avg_inc_all(1) + ms_avg_inc_all(2), ms_avg_inc_all(3) + ms_avg_inc_all(4), ...
                 ms_avg_inc_all(5), ms_avg_inc_all(6) + ms_avg_inc_all(7),  ...
                 ms_avg_inc_all(8) + ms_avg_inc_all(9), ms_avg_inc_all(10), ms_avg_inc_all(11)];
% Visualize distribution of contract valuations.
                 v_i_all_aux = repmat(delta_opt,1,NS_cons) + mu_obs;
v_i_all = permute(reshape(v_i_all_aux,[J,T,NS_cons]),[2,1,3]);             
for j = 1:J
    subplot(3,4,j)
    % Select relevant data to plot.
    v_select = reshape(squeeze(v_i_all(:,j,:)),NS_cons*T,1);
    avg = mean(v_select);
    hist(v_select);
    title_str = sprintf('Contract %d',j);
    title(title_str);
    xlim([-6;6]);
    legend({sprintf('mean = %3.2f', avg)});
end
sgtitle('Distribution of contract valuation distribution (across time and consumers)');
saveas(gcf,project_paths('OUT_FIGURES',['contract_val_dist_',file_suffix_full,'.pdf']));
% Save raw valuation data in mat-file for visualizing in Python.
% 3-dimensional matrix with time-contract-consumer structure.
save(project_paths('OUT_ANALYSIS',['contract_val_dist_',file_suffix_full,'.mat']),'v_i_all');

%% Repeat exercise with various markups on wholesale price (wholesale-plus contracts)
% Compute new price charged by consumer. 
% Markup that firm charges to consumers (in percent above wholesale spot
% price).
% Compute welfare for grid of markups.
% Because here model has to be resolved, let's focus on fewer markup levels for exploration purposes.
markup_grid_all = linspace(0,250,751)' ./100;
CS_CV_inc_all_WHP_grid = zeros(T,3,length(markup_grid_all));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alternative approach: compute markups in closed-form
% For sanity check: Compare observed retail prices to wholesale spot price in every month.
% Need to add outside good price to price in data (which is already
% normalized for outside good).
% price_retMarkup_grid = zeros(T,J,length(markup_grid));
% Compute actual change in consumer surplus for a given markup level.
CS_inc_all_WHP = zeros(T,3);
for mu_idx=1:length(markup_grid_all)
    markup = markup_grid_all(mu_idx);
    % Add markup to correctly scaled wholesale price measure.
    p_whmu = p_wh .* (1.0 + markup) .* vat_rate_wide;
    p_whmu_norm = p_whmu - repmat(p_og,1,J);
    % Compute consumer surplus by adding price differential between observed and
    % hypothetical incumbent price.
    price_change = p_whmu_norm - price_data_check;
    price_util_adjustment = (repmat(price_change,1,1,NS_cons) .* permute(repmat(price_coeff_obs,1,T,J),[2,3,1]));
    v_inc_cf_all_aux = repmat(delta_wide,[1,1,NS_cons]) + permute(mu_wide,[2,1,3]) + price_util_adjustment;
    
    v_inc_cf_all_WHP_all = squeeze(log(sum(exp(v_inc_cf_all_aux),2)));
    v_inc_cf_all_WHP_young = v_inc_cf_all_WHP_all(:,preference_shocks(1,:)'==0);
    v_inc_cf_all_WHP_old = v_inc_cf_all_WHP_all(:,preference_shocks(1,:)'==1);

    % Compute new consumer surplus for different consumer types.
    % Col 1: Averaged over all consumers
    % Col 2: Averaged only over young consumers.
    % Col 3: Averaged only over old consumers.
    CS_inc_all_WHP(:,1) = - mean(v_inc_cf_all_WHP_all ./ repmat(price_coeff_obs',T,1),2);
    CS_inc_all_WHP(:,2) = - mean(v_inc_cf_all_WHP_young ./ repmat(price_coeff_young',T,1),2);
    CS_inc_all_WHP(:,3) = - mean(v_inc_cf_all_WHP_old ./ repmat(price_coeff_old',T,1),2);
    CS_CV_inc_all_WHP = 100 * (CS_inc_all_WHP - CS_data_agg);
    CS_CV_inc_all_WHP_grid(:,:,mu_idx) = CS_CV_inc_all_WHP;   
end

% Plot size of consumer surplus change for different markup levels and consumer types.
subplot(1,3,1)
sgtitle('Change in consumer surplus for wholesale-plus contract (incumbent all)')
CS_CV_inc_all_WHP_plot = squeeze(mean(CS_CV_inc_all_WHP_grid,1));
plot(markup_grid_all, CS_CV_inc_all_WHP_plot(1,:));
title('All consumers');
ylabel('Avg change in consumer surplus per month (in EUR)')
xlabel('Markup charged on wholesale spot price')
subplot(1,3,2)
plot(markup_grid_all, CS_CV_inc_all_WHP_plot(2,:));
title('Young consumers');
ylabel('Avg change in consumer surplus per month (in EUR)')
xlabel('Markup charged on wholesale spot price')
subplot(1,3,3)
plot(markup_grid_all, CS_CV_inc_all_WHP_plot(3,:));
title('Old consumers');
ylabel('Avg change in consumer surplus per month (in EUR)')
xlabel('Markup charged on wholesale spot price')
% Find consumer surplus change closest to zero.
[CS_CV_inc_all_WHP_sol, CS_CV_inc_all_WHP_sol_idx] = min(abs(CS_CV_inc_all_WHP_plot),[],2);
CS_CV_inc_all_WHP_sol_mu = markup_grid(CS_CV_inc_all_WHP_sol_idx);
% Do actual markup comparison (of potential interest).
price_MUChange = price_retMarkup_grid(:,:,CS_CV_inc_all_WHP_sol_idx);
fprintf('Consumers are indifferent between status quo and hypothetical monopolist\n with all existing contracts when markup is %4.2f - %4.2f - %4.2f percent.\n',100*CS_CV_inc_all_WHP_sol_mu);
% Observed price time-series (in EUR) for incumbent contract.
price_cf_inc_all_WHP = [price_ret_wide(1,:)' .* 100, 100 *p_wh(:,1) .* (1+vat_rate) .* (1.0 + CS_CV_inc_all_WHP_sol_mu(1))]; 
% Save raw price time series data for formattign in Python.
% Write to table for formatting in Python.
price_labels = {'Observed','Indiff. Markup'};
cf_inc_all_id_price = table(price_cf_inc_all_WHP(:,1),price_cf_inc_all_WHP(:,2), 'VariableNames', price_labels);
writetable(cf_inc_all_id_price, project_paths('OUT_ANALYSIS',['cf_inc_all_price_id_',file_suffix_full,'.csv']));
% Exploratory graph.
close all
subplot(1,1,1)
plot(time_grid, price_cf_inc_all_WHP(:,1), time_grid, price_cf_inc_all_WHP(:,2));
title('Change in consumer surplus for wholesale-plus contract (incumbent all)')
ylabel('Monthly retail price (in EUR)')
xlabel('Time')
legend('Observed','Indifference markup')
saveas(gcf,project_paths('OUT_FIGURES',['prices_cf_inc_all_WHP_',file_suffix,'.pdf']))

% Compute hypothetical market share distribution for indifference markup.
% Compute market share predictions for case where incumbent ofers all
% contracts at observed prices.
markup_inc_all_opt = CS_CV_inc_all_WHP_sol_mu(1);
% Add markup to correctly scaled wholesale price measure.
p_whmu = p_wh .* (1.0 + markup_inc_all_opt);
p_whmu_norm = p_whmu - repmat(p_og,1,J);
% Compute consumer surplus by adding price differential between observed and
% hypothetical incumbent price.
price_change = p_whmu_norm - price_data_check;
price_util_adjustment = (repmat(price_change,1,1,NS_cons) .* permute(repmat(price_coeff_obs,1,T,J),[2,3,1]));
v_inc_cf_all_opt = repmat(delta_wide,[1,1,NS_cons]) + permute(mu_wide,[2,1,3]) + price_util_adjustment;
 
% Compute shares for indifference markup when incumbent offers all contracts. 
share_inc_all_WHP = ... numerator
    exp(v_inc_cf_all_opt) ...
    ./ ... denominator
repmat(sum(...
    exp(v_inc_cf_all_opt) ...
    ,2),[1,J,1]);

ms_avg_inc_all_WHP = squeeze(mean(mean(share_inc_all_WHP,3),1))';
ms_avg_firm_inc_all_WHP =  ... % Aggregate to firm level in order to avoid confidentiality issues.
                [ms_avg_inc_all_WHP(1) + ms_avg_inc_all_WHP(2), ms_avg_inc_all_WHP(3) + ms_avg_inc_all_WHP(4), ...
                 ms_avg_inc_all_WHP(5), ms_avg_inc_all_WHP(6) + ms_avg_inc_all_WHP(7),  ...
                 ms_avg_inc_all_WHP(8) + ms_avg_inc_all_WHP(9), ms_avg_inc_all_WHP(10), ms_avg_inc_all_WHP(11)];

%% New counterfactual #2: some consumers are put on cheapest tariff every year?
% idea:
%{
    - Think about service that puts an individual consumer on the cheapest
    contract every quarter/6 months/year?
    - Assume this program is small so that all prices remain unchanged.
    - What's the generated welfare gain? No search and switching frictions
    anymore, always on cheapest price, but potentially bad quality, etc.
    - Overall, I still expect a very large welfare gain, i.e., consumers
    should have a pretty big willingness to pay for such a service.
%}
% Construct tie series of chepeast supplier every month.
[cheapeast_price,cheapest_idx] = min(price_data_check,[],2);

% Potentially, we can modify this to only update the consumer every
% quarter?
% Potentially, discriminate between conventional and green electricity? 
% Again this would only require selecting the cheapeast among only green
% contracts.

% Create matrix for how often cheapest supplier gets adjusted.
cheap_idx_matrix = zeros(T,4);
% For simplicty, assume adjustment starts in first year of the month.
% Every month as baseline.
cheap_idx_matrix = cheapest_idx;
%% Every quarter.
cheap_idx_quarter = cheapest_idx;
% Adjust for first quarter.
cheap_idx_quarter(2) = cheap_idx_quarter(1);
for t_q=3:3:T
    cheap_quarter = cheapest_idx(t_q);
	cheap_idx_quarter(t_q) = cheap_quarter;
    cheap_idx_quarter(t_q+1:t_q+2) = cheap_quarter;
    %cheapest_idx_quarter(t_q+2) = cheap_quarter;
end
cheap_idx_matrix(:,2) = cheap_idx_quarter;
%% Every six months.
cheap_idx_halfyear = cheapest_idx;
% Adjust for first quarter.
cheap_idx_halfyear(2:5) = cheap_idx_halfyear(1);
for t_q=6:6:T
    cheap_halfyear = cheapest_idx(t_q);
	cheap_idx_halfyear(t_q) = cheap_halfyear;
    cheap_idx_halfyear(t_q+1:t_q+5) = cheap_halfyear;
    %cheapest_idx_quarter(t_q+2) = cheap_quarter;
end
cheap_idx_matrix(:,3) = cheap_idx_halfyear;
% Every year.
cheap_idx_year = cheapest_idx;
cheap_idx_year(2:11) = cheap_idx_year(1);
for t_q=12:12:T
    cheap_year = cheapest_idx(t_q);
	cheap_idx_year(t_q) = cheap_year;
    cheap_idx_year(t_q+1:min(t_q+11,T)) = cheap_year;
    %cheapest_idx_quarter(t_q+2) = cheap_quarter;
end
cheap_idx_matrix(:,4) = cheap_idx_year;

% Compute deltas and mus for cheapeast supplier path.
delta_cheap = zeros(T,NS_cons);
mu_cheap = zeros(T,NS_cons);
delta_obs_wide = reshape(log(expmvalold_backup),J,T)';
mu_obs_wide = reshape(mu_obs,[J,T,NS_cons]);

% (not ideal) but loop over months to select cheapest delta and mus.
% Here, we can adjust the frequency at which the program reoptimies
% consumers contract choice:
% 1. every month
% 2. every quarter
% 3. every 6 months
% 4. every year
for t=1:T
    cheap_idx_t = cheap_idx_matrix(t,1);
    delta_cheap(t,:) = delta_obs_wide(t,cheap_idx_t);
    mu_cheap(t,:) = squeeze(mu_obs_wide(cheap_idx_t,t,:))';
end
% Compute total utility from always being on the cheapest contract.
v_cheap_cf = delta_cheap + mu_cheap;
% Compute new consumer surplus.
CS_cheap = - mean((v_cheap_cf) ./ repmat(price_coeff_obs',T,1),2);
CS_CV_cheap = 100 * (CS_cheap - CS_data_agg);
fprintf('Average increase in consumer surplus per month and consumer (in EUR) with automatic-cheap meachanism \n  All \t Young \t Old \n %4.2f \t %4.2f \t %4.2f \n',mean(CS_CV_cheap));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF COUNTERFACTUAL COMPUTATIONS.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FORMAT WELFARE PREDICTIONS.
cf_labels = {'Reduced switching costs', 'Reduced PCW search costs','Regulator campaign','Auto. lowest price','Incumbent only conv.','Incumbent only','Incumbent only (no SC)','Incumbent only (logit corr.)','Incumbent only (logit corr. and no SC)','Incumbent all (no SC)'};
col_labels = {'Average','Non-seniors','Seniors'};

welfare_CV_comp = [ mean(CS_CV_subsidy_switch); ... # reduced switching costs
                    mean(CS_CV_perfect_info); ... # reduced PCW search cost
                    2.7 * sum(CS_CV_campaign,1); ... # for regulator campaign, we report total welfare gain for all of Flanders.
                    mean(CS_CV_cheap);... # automatic cheap mechanism
                    mean(CS_CV_inc_only_1); ... # only conventional incumbent contract is available with perfect info, and trivially no contract-specific switching costs
                    mean(CS_CV_inc_only_2); ... # only incumbent contracts are available with perfect info, but still contract-specific switching costs
                    mean(CS_CV_inc_only_2_nosc); ... only incumbent offers contracts and switching costs across contracts are also removed
                    mean(CS_CV_inc_only_2_logitJ); ... only incumbent offers contracts but logit shocks are "correcetd for"
                    mean(CS_CV_inc_only_2_logitJ_nosc); ... only incumbent offers contracts but logit shocks are "correcetd for" and all market frictions , i.e., switching costs, are eliminated.
                    mean(CS_CV_inc_all); ... # incumbent offers all available contracts without any frictions -> best-case scenario
                    ];
% Write to table for formatting in Python.
welfare_comp_table = table(welfare_CV_comp(:,1),welfare_CV_comp(:,2),welfare_CV_comp(:,3), 'RowNames',cf_labels,'VariableNames', col_labels);
writetable(welfare_comp_table, project_paths('OUT_ANALYSIS',['welfare_comp_',file_suffix_full,'.csv']),'WriteRowNames',true);

% Format output of all counterfactuals.
%% Export table with overview of observed, predicted, perfect info and switching subsidy.
counterfactual_ms = [ms_avg_firm_obs', ms_avg_firm_pred', ...
                     ms_avg_firm_subsidy', ms_avg_firm_perfect_info', ...
                     ms_avg_firm_reg_campaign', ...
                     ms_avg_firm_inc_all', ...
                     ms_avg_firm_inc_all_WHP', ...
                     ];
cf_ms_col_labels = {'Observed', 'Predicted','Switching subsidy', 'Perfect information', 'Regulator Campaign', 'Incumbent - All', 'Incumbent - All Wh+'};
row_labels = {'ECS', ...
              'EDF', ...
              'Essent', ...
              'ENINuon', ...
              'Eneco', ...
              'Lampiris', ...
              'Other'};

% Write raw data to file for table formatting in Python.
cf_ms_table = table( ...
     counterfactual_ms(:,1), counterfactual_ms(:,2), counterfactual_ms(:,3), ...
     counterfactual_ms(:,4), counterfactual_ms(:,5), counterfactual_ms(:,6), counterfactual_ms(:,7), ...
     'VariableNames', cf_ms_col_labels,'RowNames',row_labels);
writetable(cf_ms_table, project_paths('OUT_ANALYSIS',['cf_marketshares_',file_suffix_full,'.csv']));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRINT COMPARISON OF SHARES, CHURN, PCW USAGE FOR ALL COUNTERFACTUALS.

%% Check observed, predicted and counterfactual churn rates.
fprintf('\n\n :MODEL PREDICTIONS FOR CHURN RATES: \n\n');
fprintf('Average monthly churn observed in the data: \n %6.4f \n',mean(churn_rates));
fprintf('Average monthly churn predicted by model: \n %6.4f \n',mean(pred_churnAgg1d));
fprintf('Average monthly churn with switching subsidy of EUR %2.0f: \n %6.4f \n', subsidy_switch*100, mean(cr_subsidy));
fprintf('Average monthly churn when PCW search cost is cut by %2.0f percent: \n %6.4f \n', kappa_cut*100, mean(cr_perfect_info));
fprintf('Average monthly churn when regular campaign is eliminated: \n %6.4f \n',mean(cr_reg_campaign));

%% Compute effects on consumers' search behavior, i.e. usage of PCW.
pcw_use_obs = pcw_use_macro(:,1);
% Compute time series of PCW usage for observed data, predicted by model,
% switching subsidy counterfactual, and counterfactual reducing search
% costs.
PCW_use_comparison = [pred_PCWAgg1d pcw_use_subsidy_switch pcw_use_perfect_info pcw_use_reg_campaign]; % pcw_use_inc_all, pcw_use_inc_allwhp];
PCW_use_comparison = [pcw_use_obs PCW_use_comparison];
% Compute relative increase in PCW use for different scenarios in percent
PCW_change = [pcw_use_subsidy_switch pcw_use_perfect_info pcw_use_reg_campaign... pcw_use_inc_all pcw_use_inc_allwhp
    ] - repmat(pred_PCWAgg1d,1,3);
fprintf('\n\n :MODEL PREDICTIONS FOR PCW USAGE: \n\n')
fprintf('Average share of PCW users observed in the data: \n %6.4f \n',mean(PCW_use_comparison(:,1)));
fprintf('Average share of PCW users predicted by model: \n %6.4f \n',mean(PCW_use_comparison(:,2)));
fprintf('Average share of PCW users with switching subsidy of EUR %2.0f: \n %6.4f \n', subsidy_switch*100, mean(PCW_use_comparison(:,3)));
fprintf('Average share of PCW users when PCW search cost is cut by %2.0f percent: \n %6.4f \n', kappa_cut*100, mean(PCW_use_comparison(:,4)));
fprintf('Average share of PCW users when regulator campaign is eliminated: \n %6.4f \n',mean(PCW_use_comparison(:,5)));

%% Export time series of PCW usage to csv-file for graph and table formatting in Python.
% Write to table for formatting in Python.
pcw_comp_labels = {'Observed','Predicted','Switching Subsidy', 'Search Subsidy','Reg. campagign'};
cf_pcw_use_table = table(PCW_use_comparison(:,1), PCW_use_comparison(:,2), PCW_use_comparison(:,3), ...
    PCW_use_comparison(:,4), PCW_use_comparison(:,5), 'VariableNames', pcw_comp_labels);
writetable(cf_pcw_use_table, project_paths('OUT_ANALYSIS',['cf_pcw_use_',file_suffix_full,'.csv']));

pcw_change_labels = {'Switching Subsidy', 'Search Subsidy','Reg. campagign'};
cf_pcw_change_table = table(PCW_change(:,1), PCW_change(:,2), PCW_change(:,3), ...
    'VariableNames', pcw_change_labels);
writetable(cf_pcw_change_table, project_paths('OUT_ANALYSIS',['cf_pcw_change_',file_suffix_full,'.csv']));

%% Collect indifference markups in vector to output for table formatting in
% Python.
id_markup_vec = [CS_CV_inc_only_1_WHP_sol_mu'; ...
                 CS_CV_inc_only_2_WHP_sol_mu'; ...
                 CS_CV_inc_2_nosc_WHP_sol_mu'; ...
                 CS_CV_inc_2_lc_nosc_WHP_sol_mu'; ...
                 CS_CV_inc_all_WHP_sol_mu'; ...
                 repmat(mean(mean(markup_wh)),1,3)];
id_col_labels = {'id_mu_avg','id_mu_young','id_mu_senior'}
id_mu_labels = {'Inc only conv.','Inc only', 'Inc only (no SC)','Inc only (LC, no SC)','Inc all', 'Obs'};
id_mu_table = table(id_markup_vec(:,1),id_markup_vec(:,2),id_markup_vec(:,3), 'VariableNames', id_col_labels,'RowNames',id_mu_labels);
writetable(id_mu_table, project_paths('OUT_ANALYSIS',['id_mu_',file_suffix_full,'.csv']),'WriteRowNames',true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close counterfactual file and save all output.

% Save counterfactual workspace.
if efficient_GMM==0
     save(project_paths('OUT_ANALYSIS',['cf_workspace_1_',file_suffix]),'-v7.3');
elseif efficient_GMM==1
     save(project_paths('OUT_ANALYSIS',['cf_workspace_2_',file_suffix]),'-v7.3');
end
end % end loop over model types.
% Turn of log file.
diary OFF