%{

    Model file. 
    The file allows for eight different model specifications depending on how 
    psi (switching cost) and kappa (PCW search cost) are modeled.
    1. Both psi and kappa homogeneous (baseline model)
    2. psi homogeneous, kappa as a function of age. (not estimated in final version)
    3. psi as function of age, kappa homogeneous. (not estimated in final version)
    4. Both psi and kappa as function of age. (model extension)
    5. No (transactional) switching costs, psi=0. (restricted model 1, no switching costs)
    6. No PCW search costs, kappa=0. (restricted model 2, no PCW search costs)
    7. Baseline model without Hausman IV
    8. Baseline model with alternative value of logit smoother

    Inputs for the model file:
    1. theta: vector of parameters.
    2. cdid and cdindex: Nevo-style index vectors.
    3. market_shares: contract-level market share data.
    4. s_init: initial market share distribution at beginning of first period.
    5. churn_rates: aggregate churn rates.
    6. mm_data: Data for constructing micromoments.
    7. pcw_use_micro: consumer-level information on PCW usage.
    8. pcw_inst_micro: instruments to interact individual-level PCW usage prediction error with.
    9. pcw_use_macro: data about aggregate PCW usage.
    10. adv_data: data on supplier-level advertising.
    11. adv_data_mm: data on supplier-level advertising reshaped for use with micromoments.
    12. demo_data: not used in main model specification.
    13. X_reg: regressors to be used for compuation of individual-specific part of utility and in profiling of linear parameters.
    14. Z_matrix: multi-dimensional array that stacks all IVs used for the 6 different sets of moments.
    15. draws_p_normalized: Draws for simulating consumers' beliefs about prices of "unknown" products.
    16. draws_eps: Draws for logit shocks when simulating PCW usage decision.
    17. preference_shocks: consumer-level demographics and draw for green energy preference, used to compute consumer-specific part of utility.
    18. adv_shocks: simulated shocks for determining each consumer's fallback chocie set from advertising.
    19. WMatrix: GMM weighting matrix.
    20. tol_delta_c, max_iter_c: convergence criteria for contraction mapping code.

    Outputs of the model file:
    1. model_myopic: objective function value
    2.,3.,4.,5.,6.: moments matrices
        G_obs_1: BLP moments interacted with xi
        G_obs_2: Aggregate churn rate prediction error moments
        G_obs_3: Aggregate PCW usage prediction error moments
        G_obs_4: Micro-level PCW usage moments
        G_obs_5: Micro-level switching heterogeneity moments
        G_obs_6: Micro-level contract/firm choice moments (separate sets
        for four different survye respondents)
    7. predicted aggregate market shares (month-contract)
    8. predicted aggregate churn rate (time-series)
    9. predicted consumer surplus for PCW users (time-series)
    10. predicted consumer surplus for PCW non-users (time-series)
    11. predicted PCW usage rates (time-series)
    12. consumer-specific portion of utility
    13. consumer-specific price coefficients
    14. consumer-specific CCP matrices
    15. micromoments: consumer-level predictions for using PCW.
    16. micromoments: consumer-level predictions for switching ratio.
    17. consumer-specific consumer surplus

%}


function [model_myopic, ...
    G_obs_1, G_obs_2, G_obs_3, G_obs_4, G_obs_5, G_obs_6, ...
    pSharesAgg1d, pChurnAgg1d, cSurplusPCWAgg1d, cSurplusNoPCWAgg1d,pPCWAgg1d, ...
    mu_utility,p_coeff_i, ...
    ccpiPredAvg,mm_predict_active,mm_predict_sw_ratio,cSurplusI] = model_myopic(theta, cdid, cdindex, market_shares, s_init, churn_rates, ...
                                mm_data, pcw_use_micro, pcw_inst_micro, pcw_use_macro, adv_data, adv_data_mm, demo_data, ...
                                X_reg, Z_matrix, ...
                                draws_p_normalized, draws_eps, preference_shocks, adv_shocks,...
                                WMatrix, tol_delta_c, max_iter_c)

    %% Set global variables and auxiliary parameters.
    global J JF T NS_cons NS_p_eps  ...
    run_estimation SE_calc run_counterfactuals ...
    cf_subsidy_switch subsidy_switch ...
    cf_perfect_info kappa_cut ...
    cf_reg_campaign cf_inc_only_all ...
    model_type ...
    g_bar_all theta_all debug_model comp_old_sc_CCP
    % Timer for overall evaluation of the model.
    tic;
    %% Only fur debugging: Collect all parameter guesses for all iterations.
    %theta_all = [theta_all, theta];
    
    % Define model type file suffix.
    file_suffix = ['mod',num2str(model_type)];

    %% Extract parameters from vector depending on model specification.
    % These are the non-linear parameters.
    % They depend on which of the 4 model specifications we estimate.
    if model_type==1 || model_type==7 || model_type==8
        % For model with homogeneous kappa and homogeneous psi.
        % Parameters in PCW-usage model.
        theta_pcw = theta(1:3);
        % Price coefficient.
        p_coeff = theta(4);
        % Standard deviation of green coefficient.
        var_green = theta(5);
        % Interaction term: senior*incumbent
        int_age = theta(6);
        % Interaction: price*low income
        int_income = theta(7);
        % Switching cost.
        sc_coeff = theta(8);
        % Advertising coefficients.
        adv_coeff = theta(9:10);
    elseif model_type==2
        % For model with demographic-specific kappa and homogeneous psi.
        % Parameters in PCW-usage model.
        theta_pcw = theta(1:4);
        % Price coefficient.
        p_coeff = theta(5);
        % Standard deviation of green coefficient.
        var_green = theta(6);
        % Interaction term: senior*incumbent
        int_age = theta(7);
        % Interaction: price*low income
        int_income = theta(8);
        % Switching cost.
        sc_coeff = theta(9);
        % Advertising coefficients.
        adv_coeff = theta(10:11);
    elseif model_type==3
        % For model with homogeneous kappa and heterogeneous psi.
        % Parameters in PCW-usage model.
        theta_pcw = theta(1:3);
        % Price coefficient.
        p_coeff = theta(4);
        % Standard deviation of green coefficient.
        var_green = theta(5);
        % Interaction term: senior*incumbent
        int_age = theta(6);
        % Interaction: price*low income
        int_income = theta(7);
        % Switching cost.
        sc_coeff = theta(8:9);
        % Advertising coefficients.
        adv_coeff = theta(10:11);
    elseif model_type==4
        % For model with demographic-specific kappa and demographic-specific psi.
        % Parameters in PCW-usage model.
        theta_pcw = theta(1:4);
        % Price coefficient.
        p_coeff = theta(5);
        % Standard deviation of green coefficient.
        var_green = theta(6);
        % Interaction term: senior*incumbent
        int_age = theta(7);
        % Interaction: price*low income
        int_income = theta(8);
        % Switching cost.
        sc_coeff = theta(9:10);
        % Advertising coefficients.
        adv_coeff = theta(11:12);
    elseif model_type==5
        % For model with homogeneous kappa and NO psi.
        % Parameters in PCW-usage model.
        theta_pcw = theta(1:3);
        % Price coefficient.
        p_coeff = theta(4);
        % Standard deviation of green coefficient.
        var_green = theta(5);
        % Interaction term: senior*incumbent
        int_age = theta(6);
        % Interaction: price*low income
        int_income = theta(7);
        % Switching cost.
        sc_coeff = 0;
        % Advertising coefficients as constructed in model file.
        adv_coeff = theta(8:9);
    elseif model_type==6
        % For model with NO kappa and homogeneous psi.
        % Price coefficient.
        p_coeff = theta(1);
        % Standard deviation of green coefficient.
        var_green = theta(2);
        % Interaction term: senior*incumbent
        int_age = theta(3);
        % Interaction: price*low income
        int_income = theta(4);
        % Switching cost.
        sc_coeff = theta(5);
        % Advertising coefficients as constructed in model file.
        adv_coeff = theta(6:7);
    end
    
   
    
    % Compute vector of individual specific search and switching costs.
    % Compute PCW search cost for each simulated consumer.
    % Compute as sum of aggregate time-series component and individual
    % specific senior component.
    % Models with homogeneous kappa.
    if model_type==1 || model_type==3 || model_type==5 || model_type==7 || model_type==8
        if cf_reg_campaign==0
            % For estimation of model and counterfactuals where PCW search
            % cost is unchanged.
            kappa_i = exp(repmat([pcw_use_macro(:,2) pcw_use_macro(:,3) pcw_use_macro(:,5)] * theta_pcw(1:3), 1, NS_cons));
            % For information counterfactual: decrease PCW search costs uniformly
            % by kappa_cut percent.
            if cf_perfect_info==1
                kappa_i = (1.0 - kappa_cut) .* kappa_i;
            end
        elseif cf_reg_campaign==1 % For counterfactual to evaluate regulator info campaign.
            kappa_i = exp(repmat([pcw_use_macro(:,2) pcw_use_macro(:,3)] * theta_pcw(1:2), 1, NS_cons));
        end
    % Model with heterogeneous kappa.
    elseif model_type==2 || model_type==4
        if cf_reg_campaign==0
        % For estimation of the model.    
            kappa_i = exp(repmat([pcw_use_macro(:,2) pcw_use_macro(:,3) pcw_use_macro(:,5)] * theta_pcw(1:3), 1, NS_cons) ...
                             + repmat(preference_shocks(1,:),T,1) .* theta_pcw(end));
        % For information counterfactual: decrease PCW search costs uniformly
        % by kappa_cut percent.
        if cf_perfect_info==1
            kappa_i = (1.0 - kappa_cut) .* kappa_i;
        end
        elseif cf_reg_campaign==1
           % For counterfactual to evaluate regulator info campaign.
                kappa_i = exp(repmat([pcw_use_macro(:,2) pcw_use_macro(:,3)] * theta_pcw(1:2), 1, NS_cons) ...
                             + repmat(preference_shocks(1,:),T,1) .* theta_pcw(end));
        end
    elseif model_type==6    
        kappa_i = zeros(T,NS_cons);
    end
    

    %% Compute vector of heterogeneous switching costs (as function of age).
    if model_type==3 || model_type==4
        sc_coeff_vec = (sc_coeff(1) .* (1-preference_shocks(1,:)') + sc_coeff(2) .* (preference_shocks(1,:)'))';
        if cf_subsidy_switch==1
            % Take switching subsidy (absolute EUR value into account when computing
            % switching costs)
            % Transform EUR-subsidy into utils.
            sc_cf_mod = -subsidy_switch .* (p_coeff + int_income .* preference_shocks(2,:));
            % Add to actual switching cost parameter.
            sc_coeff_vec = sc_coeff_vec - sc_cf_mod;
            % Compute exponentiated switching cost to add to exponentiated contract' utility when not switching.
            sc_effective = exp(sc_coeff_vec);
            sc_effective = sc_effective .* ones(T,1);
        elseif cf_subsidy_switch==0
            % Compute exponentiated switching cost to multiply to  exponentiated contract' utility when not switching.
            sc_effective = exp(sc_coeff_vec) .* ones(T,1);
        end
    % Compute vector when switching costs are homogenous across consumers.
    elseif model_type==1 || model_type==2 || model_type==5 || model_type==6 || model_type==7 || model_type==8
        if cf_subsidy_switch==1
            % Take switching subsidy (absolute EUR value into account when computing
            % switching costs)
            % Transform EUR-subsidy into utils.
            sc_cf_mod = -subsidy_switch .* (p_coeff + int_income .* preference_shocks(2,:));
            % Add to actual switching cost parameter.
            sc_coeff = sc_coeff - sc_cf_mod;
            % Compute exponentiated switching cost to add to exponentiated contract' utility when not switching.
            sc_effective = exp(sc_coeff);
            sc_effective = sc_effective .* ones(T,1);
        elseif cf_subsidy_switch==0
            % Compute exponentiated switching cost to add to exponentiated contract' utility when not switching.
            sc_effective = exp(sc_coeff);
            sc_effective = sc_effective .* ones(T,NS_cons);
        end
    end
    % Reshape search and switching costs into long-format.
    % PCW search costs needs to be exponentiated for C code.
    % Order: months - simulated consumers.
    exp_kappa_i_long = exp(reshape(kappa_i',NS_cons*T,1));
    % Switching costs are already exponentiated above.
    % Order (shouldn't matter for most models, since switching cost are
    % constant over time): months - simulated consumers.
    exp_psi_i_long = reshape(sc_effective',NS_cons*T,1);
    
    % When simulating one incumbent that offers all contracts without frictions, set
    % both psi and kappa to zero (or exponentiated versions equal to 1).
    if cf_inc_only_all==1
        exp_kappa_i_long = ones(NS_cons*T,1);
        exp_psi_i_long = ones(NS_cons*T,1);
    end
    % Compute individual-specific price coefficient using
    % individual-specific income deviation.
    p_coeff_i = repmat(p_coeff,NS_cons,1) + int_income .* preference_shocks(2,:)';
    % END OF CONSTRUCTING INDIVIDUAL-SPECIFIC PARAMETER VECTORS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Compute exogenous awareness sets through advertising.
    % Small version with just two parameters for debugging.
    adv_index = adv_coeff(1) .* adv_data(:,:,:,1) ...
                + adv_coeff(2) .* adv_data(:,:,:,2);
    % Predicted probability of a consumer being aware.
    prob_aware_probit = normcdf(adv_index);
    % Main specification has probit for this process.
    prob_aware = prob_aware_probit;
   
    % Compare to awareness shocks to simulated choice sets.
    % Simulate whether in a given months (rows), a consumer (height), is
    % aware of a firm (columns).
    aware_sim = (prob_aware>1.0*adv_shocks);
    % For debugging: What happens when everybody is fully informed.
    % aware_sim = ones([T,7,NS_cons]);
    
    % Expand awareness from firm to contract level: replicate columns for
    % firm depending on how many contracts each firm has, i.e., we go from
    % 7 to 11 columns.
    info_index = cat(2,repmat(aware_sim(:,1,:,:),1,2,1,1), repmat(aware_sim(:,2,:,:),1,2,1,1), repmat(aware_sim(:,3,:,:),1,1,1,1), repmat(aware_sim(:,4,:,:),1,2,1,1), repmat(aware_sim(:,5,:,:),1,2,1,1), repmat(aware_sim(:,6,:,:),1,1,1,1),repmat(aware_sim(:,7,:,:),1,1,1,1)); 
    % Expand this into CCP-matrix compatible format. 
    % Dim 1: previous contract
    % Dim 2: current contract choice
    % Dim 3: month
    % Dim 4: simulated consumer
    info_index = permute(repmat(permute(info_index,[2,5,1,3,4]),[1,J,1,1,1]),[2,1,3,4,5]);
    % Adjust info matrix for existing firm always being available.
    % Construct index matrix for existing firms' contracts.
    ec_index = [repmat([ones(2,1);zeros(9,1)],1,2), ...
                repmat([zeros(2,1);ones(2,1);zeros(7,1)],1,2), ...
                [zeros(4,1);ones(1,1);zeros(6,1)], ...
                repmat([zeros(5,1);ones(2,1);zeros(4,1)],1,2), ...
                repmat([zeros(7,1);ones(2,1);zeros(2,1)],1,2), ...
                [zeros(9,1);ones(1,1);zeros(1,1)], ...
                [zeros(10,1);ones(1,1)]           
                ];
    avg_contracts_through_adv = mean(mean(sum(aware_sim,2)));
  
    % Adjust info index matrix for being always aware of existing contract.
    % Dimension assignment as for info_index object above.
    info_index(repmat(ec_index,[1,1,T,NS_cons])==1) = 1;
    
    % Reshape matrix into long format according to following order:
    % Month - consumer - lagged choice - current choice
    info_index_c = permute(info_index,[2,1,4,3]);
    % If memory becomes an issue, this could be cast as integer.
    % Since it is not passed around much, doube should be fine in terms of
    % speed.
    info_index_long = double(reshape(info_index_c,J*J*NS_cons*T,1));
      
    % END OF AWARENESS-THROUGH-ADVERTISING-SIMULATION.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Compute mu and exp(mu) terms, ie, consumer specific utility term.
    % New version where green valuation is log-normally distributed.
    % Linear parameter captures lower end of log-normal distribution.
    % mu-term only contains positive deviation from lower bound.
    % First compute mu-part associated with demographics.
    mu_utility_demo =  [X_reg(:,1), X_reg(:,2)] * (repmat([int_age; int_income],1,NS_cons) .* preference_shocks(1:2,:)); 
    mu_utility_green = exp(X_reg(:,3) * (repmat(var_green,1,NS_cons) .* preference_shocks(3,:))); 
    mu_utility_green(mu_utility_green==1) = 0;
    mu_utility = mu_utility_demo + mu_utility_green;
    % Old Version where green valuation is normally distributed.
    % mu_utility =  [X_reg(:,2), X_reg(:,3), X_reg(:,4)] * (repmat([int_age; int_income; var_green],1,NS_cons) .* preference_shocks(1:3,:)); 
    expmu = exp(mu_utility);
    
    % Multiply price belief difference with price coefficient.
    % Expand consumer-specific price vector to conform with price belief draws.
    p_coeff_expand = repmat(kron(p_coeff_i,ones(NS_p_eps*J,1)),T,1);
    % Remember that price belief draws are already exponentiated:
    draws_pbelief_alpha = draws_p_normalized .^ p_coeff_expand;
    
    % Load old exponentiated mean utilities in vector format from previous iteration.
    % This includes mean utility of outside good normalized to 0
    % explicitly.
    if cf_inc_only_all==0
        load(project_paths('OUT_ANALYSIS',['expmvalold_',file_suffix,'.mat']),'expmvalold');
    elseif cf_inc_only_all==1
        load(project_paths('OUT_ANALYSIS',['expmvalold_cfincall_aux_',file_suffix,'.mat']),'expmvalold');
    end
    % END OF PREPARING MODEL INGREDIENTS FOR CONTRACTION MAPPING.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    if run_estimation==1 || SE_calc==1
        contrMap_dummy = 1;
    elseif run_estimation==0 || run_counterfactuals==1
        contrMap_dummy = 0;
    end
    
    %% Call contraction mapping code in C/Mex.
    
    % Timer for contraction mapping only.
    tic;
        
    % Safety measure on initial conditions.
    % If some of the demographic-type-product specific initial market shares are zero,
    % set the zeros to a tiny number. Depending on how consumer types are defined,
    % there can be a small number of zeros in the initial conditions. 
    % We experimented with different numbers here and it doesn't lead to significant changes.
    s_init_corr = reshape(s_init,J, NS_cons);
    mshares_wide = reshape(market_shares,J,T);
    s_init_agg = mean(s_init_corr,2);
    s_init_corr(s_init_corr==0) = 0.005;
    s_init_corr = s_init_corr ./ sum(s_init_corr,1);
    s_init_agg = mean(s_init_corr,2);
    s_init = reshape(s_init_corr,J*NS_cons,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if comp_old_sc_CCP==1 % only for computing CCP matrices for old SC for referee replies.
        [expmvalold, pPCWUsagei, pSharesi,pSharesiPCW, pSharesiNoPCW, pChurni, cSurplusPCWi, cSurplusNoPCWi, ccpiPred, mktnorm, mktiter] = ...
            contrMap_ES_old_SC( ...
            tol_delta_c, max_iter_c, ...
            expmvalold, expmu, cdindex, market_shares, s_init, ...
            info_index_long, draws_pbelief_alpha, draws_eps, ...
            exp_psi_i_long, exp_kappa_i_long,contrMap_dummy);
    else % this is what we use in the final version of the paper.
        if model_type==8  % for a larger value of logit smoother (only used in referee response)  
            [expmvalold, pPCWUsagei, pSharesi,pSharesiPCW, pSharesiNoPCW, pChurni, cSurplusPCWi,  cSurplusNoPCWi, ccpiPred, mktnorm, mktiter] = ...
                contrMap_ES_alt_LS( ...
                tol_delta_c, max_iter_c, ...
                expmvalold, expmu, cdindex, market_shares, s_init, ...
                info_index_long, draws_pbelief_alpha, draws_eps, ...
                exp_psi_i_long, exp_kappa_i_long,contrMap_dummy);
        else % used for final version of the paper.
                [expmvalold, pPCWUsagei, pSharesi,pSharesiPCW, pSharesiNoPCW, pChurni, cSurplusPCWi,  cSurplusNoPCWi, ccpiPred, mktnorm, mktiter] = ...
                contrMap_ES( ...
                tol_delta_c, max_iter_c, ...
                expmvalold, expmu, cdindex, market_shares, s_init, ...
                info_index_long, draws_pbelief_alpha, draws_eps, ...
                exp_psi_i_long, exp_kappa_i_long,contrMap_dummy);
        end
    end
    cm_time = toc;
    
    if run_estimation==1
        % Save updated mean utilities for loading in next iteration.
        save(project_paths('OUT_ANALYSIS',['expmvalold_',file_suffix,'.mat']),'expmvalold');
    end
    
    % At this point, we have mean utilities, and predicted shares
    % and churn rates for each individual. In addition we have predictions for PCW usage and consumer surplus. 

% Safety check for NaNs in expmvalold.
checkNaNexpmvalold =sum(isnan(expmvalold));
if checkNaNexpmvalold>0
    fprintf('!!WARNING!! NaNs in expmvalold. How many? %d.\n',checkNaNexpmvalold);
    % Load mean values from logit model if contraction mapping results in NaNs for deltas.
    % We never encountered this for this model.
    load(project_paths('OUT_ANALYSIS',['expmvalold_logit.mat']));
    save(project_paths('OUT_ANALYSIS',['expmvalold_',file_suffix,'.mat']),'expmvalold');
    % If expmvalold contain NaNs, exit function prematurely and reset delta values.
    model_myopic = 1E12;
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reshape output from C code.

% 1. Reshape market share predictions.
% Individual specific shares.
% 3D-reshaping.
pShares3d = reshape(pSharesi,NS_cons*J,T); % products x consumers - months
pShares3d = reshape(pShares3d,[J,NS_cons,T]); % products - consumers - months
pShares3d = permute(pShares3d,[3,1,2]);
% 2D-reshaping.
pShares2d = reshape(permute(pShares3d,[2,1,3]),T*J,NS_cons);

% 3D-reshaping.
pSharesPCW3d = reshape(pSharesiPCW,NS_cons*J,T); % products x consumers - months
pSharesPCW3d = reshape(pSharesPCW3d,[J,NS_cons,T]); % products - consumers - months
pSharesPCW3d = permute(pSharesPCW3d,[3,1,2]);
% 2D-reshaping.
pSharesPCW2d = reshape(permute(pSharesPCW3d,[2,1,3]),T*J,NS_cons);

% 3D-reshaping.
pSharesNoPCW3d = reshape(pSharesiNoPCW,NS_cons*J,T); % products x consumers - months
pSharesNoPCW3d = reshape(pSharesNoPCW3d,[J,NS_cons,T]); % products - consumers - months
pSharesNoPCW3d = permute(pSharesNoPCW3d,[3,1,2]);
% 2D-reshaping.
pSharesNoPCW2d = reshape(permute(pSharesNoPCW3d,[2,1,3]),T*J,NS_cons);

% Safety check for NaNs in shares conditional on PCW usage.
% We never encountered this issue for our final model specifications.
checkNaNSharesPCW =sum(sum(sum(isnan(pSharesPCW3d))));
checkNaNSharesNoPCW =sum(sum(sum(isnan(pSharesNoPCW3d))));
if checkNaNSharesPCW>0
    fprintf('!!WARNING!! NaNs in pSharesiPCW. How many? %d.\n',checkNaNSharesPCW);
    % If share predictions contain NaNs, exit function prematurely.
    model_myopic = 1E12;
    return;
end
if checkNaNSharesNoPCW>0
    fprintf('!!WARNING!! NaNs in pSharesiNoPCW. How many? %d.\n',checkNaNSharesNoPCW);
    % If share predictions contain NaNs, exit function prematurely.
    model_myopic = 1E12;
    return;
end

% Aggregate shares (matched to observed shares).
pSharesAgg2d = mean(pShares3d,3);
pSharesAgg1d = mean(pShares2d,2);
% Check matching of aggregate market shares.
cCheckShares = abs(pSharesAgg1d - market_shares);

% Construct matrix of lagged individual specific shares.
s_init2d = reshape(s_init, J,NS_cons);
s_init3d = reshape(s_init,[1,J,NS_cons]);
pSharesLag2d = [s_init2d;pShares2d(1:(T-1)*J,:)];
pSharesLag3d = cat(1,s_init3d,pShares3d(1:end-1,:,:));

% 2. Reshape churn rates.
% Individual specific churn rates.
% 3D-reshaping.
pChurn3d = reshape(pChurni,NS_cons*J,T); % products x consumers - months
pChurn3d = reshape(pChurn3d,[J,NS_cons,T]); % products - consumers - months
pChurn3d = permute(pChurn3d,[3,1,2]);
% Correct NaNs in churn rate predictions in first period for contracts 8
% and 9. This is a work-around because some firms do not have every consumer type in the
% initial period, so consumer level churn rates are not defined in period 1; therefore, replace
% with 0.
pChurn3d(1,8,isnan(pChurn3d(1,8,:))) = 0;
pChurn3d(1,9,isnan(pChurn3d(1,9,:))) = 0;
% 2D-reshaping.
pChurn2d = reshape(permute(pChurn3d,[2,1,3]),T*J,NS_cons);

% Aggregate churn time series (matched to observed churn rates).
% Weight product-specific churn with lagged product share.
pChurnAgg2d = squeeze(sum(pChurn3d .* pSharesLag3d,2));
pChurnAgg1d = mean(pChurnAgg2d,2);
% Check matching of aggregate market Churn.
cCheckChurn = abs(pChurnAgg1d - churn_rates);

% 3. Reshape PCW usage.
% Individual specific PCW usage.
% 3D-reshaping.
pPCW3d = reshape(pPCWUsagei,NS_cons*J,T); % products x consumers - months
pPCW3d = reshape(pPCW3d,[J,NS_cons,T]); % products - consumers - months
pPCW3d = permute(pPCW3d,[3,1,2]);
% 2D-reshaping.
pPCW2d = reshape(permute(pPCW3d,[2,1,3]),T*J,NS_cons);

% Aggregate PCW time series (matched to observed PCW rates).
% Weight product-specific PCW with lagged product share.
pPCWAgg2d = squeeze(sum(pPCW3d .* pSharesLag3d,2));
pPCWAgg1d = mean(pPCWAgg2d,2);
% Check matching of aggregate market PCW.
cCheckPCW = pPCWAgg1d - pcw_use_macro(:,1);

% 4. Reshape consumer surplus
% Conditional on PCW Usage. 
% 3D-reshaping.
cSurplusPCW3d = reshape(cSurplusPCWi,NS_cons*J,T); % products x consumers - months
cSurplusPCW3d = reshape(cSurplusPCW3d,[J,NS_cons,T]); % products - consumers - months
cSurplusPCW3d = permute(cSurplusPCW3d,[3,1,2]);
% 2D-reshaping.
cSurplusPCW2d = reshape(permute(cSurplusPCW3d,[2,1,3]),T*J,NS_cons);
% Aggregate PCW time series (matched to observed PCW rates).
% Weight product-specific PCW with lagged product share.
% Scale consumer surplus in unit of price coefficient (here: 100 EUR).
cSurplusPCWAgg2d = -squeeze(sum(cSurplusPCW3d .* pSharesLag3d,2)) ./ repmat(p_coeff_i',T,1);
cSurplusPCWAgg1d = mean(cSurplusPCWAgg2d,2);

% Conditional on no PCW Usage. 
% 3D-reshaping.
cSurplusNoPCW3d = reshape(cSurplusNoPCWi,NS_cons*J,T); % products x consumers - months
cSurplusNoPCW3d = reshape(cSurplusNoPCW3d,[J,NS_cons,T]); % products - consumers - months
cSurplusNoPCW3d = permute(cSurplusNoPCW3d,[3,1,2]);
% 2D-reshaping.
cSurplusNoPCW2d = reshape(permute(cSurplusNoPCW3d,[2,1,3]),T*J,NS_cons);

% Compute consumer specific surplus for each period and simulated consumer.
cSurplusI = -squeeze(sum( ...
    (cSurplusPCW3d .* pPCW3d + cSurplusNoPCW3d .* (1.0-pPCW3d))...
    .* pSharesLag3d,2)) ./ repmat(p_coeff_i',T,1);


% Aggregate PCW time series (matched to observed PCW rates).
% Weight product-specific PCW with lagged product share.
cSurplusNoPCWAgg2d = -squeeze(sum(cSurplusNoPCW3d .* pSharesLag3d,2)) ./ repmat(p_coeff_i',T,1);
cSurplusNoPCWAgg1d = mean(cSurplusNoPCWAgg2d,2);

% Weighted average of surplus (across PCW users and non-users).
% This can be reconsidered, how to exactly do the averaging of consumer
% surplus, once we have decided on the new coutnerfactuals to run.
cSurplusTotalAgg = cSurplusNoPCWAgg1d .* (1.0-pPCWAgg1d) + pPCWAgg1d .* cSurplusPCWAgg1d;

% Check matching of aggregate market PCW.
cSurplusDiff = cSurplusPCWAgg2d - cSurplusNoPCWAgg2d;

% Reshape individual-specific CCP matrices for each market.
% This reshapes in the following way:
% Rows: Lagged choice
% Columns: Current choice.
% 3th dimension: consumer types.
% 4th dimension: markets/months.
ccpiPred4d = permute(reshape(ccpiPred,[J,J,NS_cons,T]),[2,1,3,4]);
% Average over months.
% Here, it doesn't make much sense to average over consumer types
% Instead we  present exemplary consumer types' CCP matrices averaged over time.
ccpiPredAvg = mean(ccpiPred4d,4);

%% END RESHAPING CODE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if run_estimation==1
%     fprintf('Average NUMBER of firms that consumers\nare aware of through advertising: %.2f. \n', avg_contracts_through_adv);
%     fprintf('Average share of PCW users: %.4f \n', mean(pPCWAgg1d));
%     fprintf('Predicted average churn rate: %.4f \n', mean(pChurnAgg1d));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Unpack instruments array.
Z_1 = Z_matrix{1}; % BLP moments
Z_2 = Z_matrix{2}; % aggregate churn rate moments
Z_3 = Z_matrix{3}; % aggregate PCW usage moments
Z_4 = Z_matrix{4}; % micro-level PCW usage
Z_5 = Z_matrix{5}; % individual level switching heterogeneity
Z_61 = Z_matrix{6}; % Micro-level contract choice, group 1
Z_62 = Z_matrix{7}; % Micro-level contract choice, group 2
Z_63 = Z_matrix{8}; % Micro-level contract choice, group 3
Z_64 = Z_matrix{9}; % Micro-level contract choice, group 4
% Extract number of moment conditions in each moment group.
n_mc_1 = size(Z_1,2);
n_mc_2 = size(Z_2,2);
n_mc_3 = size(Z_3,2);
n_mc_4 = size(Z_4,2);
n_mc_5 = size(Z_5,2);
% Careful here: the following moment conditions are computed for each of
% the 11 (or 7) choice prediction errors.
n_mc_61 = size(Z_61,2); % * J because prediction error on contract level.
n_mc_62 = size(Z_62,2); % * 7 because prediction error on firm level.
n_mc_63 = size(Z_63,2); % * J because prediction error on contract level.
n_mc_64 = size(Z_64,2); % * 7 because prediction error on firm level.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute linear parameters
% Only first set of moments (classical BLP moments) are function of
% linear parameters, therefore, all other moments drop out from this
% profiling.
% Regress delta on exogenous linear regressors.
% This is different than in BLP for several reasons:
% Need only subset of full weighting matrix.
% -> extract relevant portion.
WM1 = WMatrix(1:n_mc_1,1:n_mc_1);
% - Our delta vector contains normalized mean utilitiy of outside good:
% -> drop outside good rows before running regression.
% - Our mean price coefficient is a nonlinear parameter:
% -> subtract mean price effect before running regression.
% -> drop price column from regressors.
delta_net_price = log(expmvalold) - p_coeff .* X_reg(:,2);
% Drop price data since price coefficient is now explicitly in contraction mapping.
X1_lin = X_reg(:,[1,3:end]);
% Drop outside good regressor and instruments.
og_idx = repmat([zeros(10,1);1],T,1);
X1_lin(og_idx==1,:) = [];
Z1_ig = Z_1(logical(og_idx==0),:);
delta_net_price(og_idx==1,:) = [];
% Update linear parameter vector only during estimation.
if run_estimation==1 % for estimation.
    P_z = Z1_ig*WM1*Z1_ig';
    theta_linear = (X1_lin'*P_z*X1_lin)\(X1_lin'*P_z*delta_net_price);
    save(project_paths('OUT_ANALYSIS',['thetalinear_',file_suffix,'.mat']),'theta_linear');
elseif SE_calc==1 % for SE computation.
    load(project_paths('OUT_ANALYSIS',['thetalinear_SE_aux_',file_suffix,'.mat']),'theta_linear');
else % for counterfactuals and other prediction exercises.
    load(project_paths('OUT_ANALYSIS',['thetalinear_',file_suffix,'.mat']),'theta_linear');
end
% Predict BLP-xi.
xi = delta_net_price - X1_lin*theta_linear; 
% END OF PROFILING OUT OF LINEAR PARAMETERS.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Only for debugging and model diagnostics. 
% This does not affect the estimation.
if debug_model==1
    % Extract the relevant parameters.
    theta_price = [p_coeff; int_income];
    theta_incumbent = [theta_linear(1); theta_linear(1) + int_age];
    theta_green = [theta_linear(2); abs(var_green)];
    % Magnitude of PCW usage costs.
    kappa_hat_young = exp([pcw_use_macro(:,2) pcw_use_macro(:,3) pcw_use_macro(:,5)] * theta_pcw(1:3)); 
    if model_type==2 || model_type==4
        kappa_hat_senior = exp([pcw_use_macro(:,2) pcw_use_macro(:,3) pcw_use_macro(:,5)] * theta_pcw(1:3) + theta_pcw(4)); 
    elseif model_type==1 || model_type==3 || model_type==5 || model_type==7 || model_type==8
        kappa_hat_senior = kappa_hat_young;
    end
    fprintf('Estimates for price coefficients for mean price and income-interaction coefficient: \n %6.4f \n %6.4f \n %6.4f \n', ...
            theta_price(1),theta_price(2))    
    
    % Compute monetary WTP for different contract characteristics.
    kappa_hat_young_EUR = -100 * kappa_hat_young ./ p_coeff(1);
    kappa_hat_senior_EUR = -100 * kappa_hat_senior ./ p_coeff(1);
    % Compute average estimate of PCW search cost.
    kappa_hat_avg_EUR = kappa_hat_young_EUR .* (1.0-0.189) +  kappa_hat_senior_EUR .* 0.189;
    kappa_out_EUR = [kappa_hat_avg_EUR, kappa_hat_young_EUR, kappa_hat_senior_EUR];
    
    % Incumbent preference in EUR: young and senior consumers.
    inc_EUR = -100 .* theta_incumbent ./ repmat(theta_price(1),2,1);
     % Green preferences in EUR.
    % Compute mean and variance of consumer-specific upward preference deviation.
    green_dev_mean = exp(0.5 * theta_green(2)^2);
    green_dev_var = exp( 2* theta_green(2)^2 - 1);
    % With lognormal distribution.
    green_EUR_lognorm = -100 ./ p_coeff .* (theta_green(1) + green_dev_mean);
    green_lower_EUR_lognorm = -100 ./ theta_price(1) .* theta_green(1);
    % With normal distribution.
    green_EUR_norm = -100 .* theta_green(1) ./ p_coeff;
    % Switching costs in EUR.
    sc_EUR = -100 * sc_coeff ./ p_coeff;
       
    %% Print implied monetary WTP for various parameters.
    % Awareness through advertising.
    fprintf('Average number of firms that consumers\nare aware of through advertising: %.2f. \n', avg_contracts_through_adv);
    % Print mean statistics on PCW usage costs.
    fprintf('Mean (over time) PCW usage costs for average, young and senior consumer, respectively:\n%4.2f EUR\n%4.2f EUR\n%4.2f EUR\n',mean(kappa_out_EUR(:,1)), mean(kappa_out_EUR(:,2)),mean(kappa_out_EUR(:,3)));
    if model_type==1 || model_type==2 || model_type==7 || model_type==8
        fprintf('Implied magnitude of switching costs for mean income consumer: %4.0f EUR.\n',...
                sc_EUR(1));
    elseif model_type==3 || model_type==4
        fprintf('Implied magnitude of switching costs for mean income young and senior consumer, respectively: %4.0f EUR and %4.0f EUR.\n',...
                sc_EUR(1), sc_EUR(2));
    end
    fprintf('Implied magnitude of incumbent preference for non-seniors (mean income consumer) : %4.2f EUR \n', ...
            inc_EUR(1));
    fprintf('Implied magnitude of incumbent preference for seniors (mean income consumer) :%4.2f EUR \n', ...
            inc_EUR(2));
    fprintf('Implied magnitude of mean green preference for mean income consumer (with normal distribution): %4.2f EUR \n', ...
            green_EUR_norm(1));
    fprintf('Implied magnitude of mean green preference for mean income consumer (with lognormal distribution) and lower bound: %4.2f EUR \n %4.2f EUR \n', ...
            green_EUR_lognorm(1), green_lower_EUR_lognorm);
    
    % Graphs for churn rates and PCW usage.
    time_grid = linspace(1,T,T)';
    subplot(2,2,1)
    set(gcf, 'PaperSize', [20 20]);
    plot(time_grid,pChurnAgg1d,time_grid,churn_rates);
    legend('Predicted churn rates','Observed churn rates');
    title('Churn Rates GOF (over time)');
    % Sanity check on predicted PCW usage.
    subplot(2,2,2)
    plot(time_grid,pPCWAgg1d,time_grid,pcw_use_macro(:,1));
    legend('Predicted PCW usage','Observed PCW usage');
    title('PCW Usage GOF (over time)');
    %saveas(gcf,project_paths('OUT_FIGURES',['gof_debug_',file_suffix,'.pdf']));
    
    %% Analyze distribution of green electricity coefficient.
    gc_grid = linspace(-5,5,100)';
    X_EUR = 5;
    coeff_X_EUR = - X_EUR .* theta_price(1) ./ 100;
    % Compute distribution of green electricity coefficient.
    % Version where green preference is modeled as normal distribution.
    gc_dist_norm = normpdf(gc_grid, theta_green(1),theta_green(2));
    green_EUR_dist_norm = -100 .* gc_grid ./ theta_price(1);
    subplot(2,2,3);
    green_wtp_norm = plot(green_EUR_dist_norm,gc_dist_norm);
    title('WTP for green electricity (normal)')
    xlabel('WTP (in EUR per month)');
    ylabel('Density');
    % saveas(green_wtp_1,project_paths('OUT_FIGURES',['wtp_green_1_',file_suffix,'.pdf']));
    sgp_norm = 1.0 - normcdf(0, theta_green(1),theta_green(2));
    fprintf('Consumer share with positive WTP for green electricity (under normal distribution): %6.4f.\n',sgp_norm);
    sgl_X_norm = 1.0 - normcdf(coeff_X_EUR, theta_green(1),theta_green(2));
    fprintf('Consumer share with WTP for green electricity larger than %d EUR (normal distribution): %6.4f.\n',X_EUR, sgl_X_norm);
    % These numbers look quite large, not sure we want to mention this in the
    % discussion of results.
    % Currently, this number looks very large to me and is another manifestation
    % of our current bad model fit, so I wouldn't stress it too much.
    fprintf('Share of consumers with positive willingness to pay for green electricity:\n %.4f\n', sgp_norm);
    fprintf('Share of consumers with willingness to pay for green electricity of more than 5 EUR:\n %.4f\n', sgl_X_norm);

    % Compute distribution of green electricity coefficient.
    % Version where green preference is modeled as lognormal distribution.
    gc_dist_logn = lognpdf(gc_grid, 0,sqrt(green_dev_var));
    green_EUR_dist_logn = -100 .* gc_grid ./ theta_price(1) + green_lower_EUR_lognorm;
    subplot(2,2,4) 
    %set(gcf, 'PaperSize', [12 24]);
    green_wtp_logn = plot(green_EUR_dist_logn,gc_dist_logn);
    title('WTP for green electricity (lognormal)')
    xlabel('WTP (in EUR per month)');
    ylabel('Density');
    %saveas(gcf,project_paths('OUT_FIGURES',['wtp_green_debug_',file_suffix,'.pdf']));
    % Save figure with model diagnostics.
    saveas(gcf,project_paths('OUT_FIGURES',['model_debug_',file_suffix,'.pdf']));

    
    sgl_X_logn = 1.0 - logncdf(coeff_X_EUR-theta_green(1), 0,theta_green(2));
    sgp_logn = 1.0 - logncdf(0-theta_green(1), 0,theta_green(2));
    fprintf('Share of consumers with positive willingness to pay for green electricity (lognormal dist):\n %.4f\n', sgp_logn);
    fprintf('Share of consumers with willingness to pay for green electricity of more than %d EUR (lognormal dist):\n %.4f\n', X_EUR, sgl_X_logn);
    
    % Average goodness of fit statistics for PCW usage and churn rates.
    fprintf('Average share of PCW users (observed): %.4f \n', mean(pcw_use_macro(:,1)));
    fprintf('Average churn rate (observed): %.4f \n', mean(churn_rates));
    fprintf('Average share of PCW users (predicted): %.4f \n', mean(pPCWAgg1d));
    fprintf('Average churn rate (predicted): %.4f \n\n', mean(pChurnAgg1d));
    fprintf('Average prediction error for PCW usage rate: %.4f \n', mean(pcw_use_macro(:,1)-pPCWAgg1d));
    fprintf('Average prediction error for churn rate: %.4f \n', mean(pChurnAgg1d));
    
    fprintf('Median share of PCW users (observed): %.4f \n', median(pcw_use_macro(:,1)));
    fprintf('Median churn rate (observed): %.4f \n', median(churn_rates));
    fprintf('Median share of PCW users (predicted): %.4f \n', median(pPCWAgg1d));
    fprintf('Median churn rate (predicted): %.4f \n\n', median(pChurnAgg1d));
    fprintf('Median prediction error for PCW usage rate: %.4f \n', median(pcw_use_macro(:,1)-pPCWAgg1d));
    fprintf('Median prediction error for churn rate: %.4f \n', median(pChurnAgg1d));
    
    fprintf('Average RMSE for PCW usage rate: %.4f \n', sqrt(mean((pcw_use_macro(:,1)-pPCWAgg1d).^2)));
    fprintf('Average RMSE for churn rate: %.4f \n', sqrt(mean((churn_rates-pChurnAgg1d).^2)));
    fprintf('\n\n!!End of diagnostics!!\n\n');
    % Look at wide format of implied expmval matrix.
    expmval_wide =reshape(expmvalold,J,T)';
    incval = sum(expmval_wide,2);
end % end of debugging and diagnostics model code. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    
    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start construction of moments.
% Order of moments is important since it needs to be consistent with the
% weighting matrix format.
% Block 1: Classical BLP moments xi * Z
% Block 2: Match aggregate churn rate moments
% Block 2a: Match heterogeneity in switching moments (based on micro data)
% Block 3: Match aggregate PCW usage moments.
% Block 4: Match micro-level PCW usage moments
% Block 5: Match micro-level contract/firm choices for the four different
% survey groups.
% Block 6: Same as 5 but interacted with advertising levels.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Block 1 of moment conditions: xi * Z
G_obs_1 = repmat(xi,1,n_mc_1) .* Z1_ig;
% end block 1 of moment conditions.

%% Block 2 of moment conditions: Aggregate churn rates.
% Compute churn rate prediction errors for the GMM objective function.
G_obs_2 = repmat((pChurnAgg1d - churn_rates),1,n_mc_2) .* Z_2;
% end block 2a of moment conditions.
% Construct PCW usage macro moments (analogous to churn rate moments).
pcw_error_macro = pPCWAgg1d - pcw_use_macro(:,1);
% Interact macro PCW-error with instruments.
G_obs_3 = repmat(pcw_error_macro,1,n_mc_3) .* Z_3;
% end block 2b of moment conditions.

%% Block 3: Micromoments based on survey responses.
% Start with contract and firm choice moments on individual level.
% (1) Extract observed contract choices from cell constructed in main file.
mm_full_known = mm_data{1};
mm_full_unknown =  mm_data{2}; 
mm_partial_known = mm_data{3}; 
mm_partial_unknown = mm_data{4}; 
% (2) Compute predicted choice probabilities for each consumer type and
% contract/firm.
% Indicators for years of our sample.
time_index = [ones(11,1); 2.*ones(12,1); 3.*ones(12,1); 4.*ones(12,1); 5.*ones(6,1)];
% Demographics index.
% Create income bins to match simulated consumers to survey data groups.
% IMPORTANT: simulated income needs to be measured in 1000-EUR!
inc_type_aux = (preference_shocks(4,:) < 1.5) ...
                + 2.* (preference_shocks(4,:) >= 1.5 & preference_shocks(4,:) < 2.5 ) ...
                + 3.* (preference_shocks(4,:) >= 2.5 &preference_shocks(4,:) < 3.8) ...
                + 4.* (preference_shocks(4,:) >= 3.8); 
% Construct composite demographics index consisting of senior-income
% characteristics.
dd_i = (4.* preference_shocks(1,:) + inc_type_aux )';
    
% (a) For fully informed consumers with contract type observed.
% These are individual contract choice probabilities for the simulated
% consumers. (months in rows, contracts in columns, consumers in height)
si_fa_mm = cat(2,repmat(time_index,[1,1,NS_cons]), pSharesPCW3d);
mm_predict_fa_contract = zeros(5*8,J);

% Average separately for each consumer type.
% Average over months within each year.
% Check whether we should use only September month for these predictions?
for c=1:8
    % Extract predicted market shares for each year.
    ms_mm_aux_1 = si_fa_mm(1:11,2:end,dd_i==c); 
    ms_mm_aux_2 = si_fa_mm(12:23,2:end,dd_i==c); 
    ms_mm_aux_3 = si_fa_mm(24:35,2:end,dd_i==c); 
    ms_mm_aux_4 = si_fa_mm(36:47,2:end,dd_i==c); 
    ms_mm_aux_5 = si_fa_mm(48:end,2:end,dd_i==c); 
    % These weights ensure that some consumers may be more likely to use PCW,
    % therefore, we weight each consumer with its probability of using the PCW.
    % Make sure these sum to one and that the resulting mm_predict objects sum to one across contracts.
    % Put differently, these weights indicate how important a specific
    % consumer is within the group of fully informed consumers.
    weights_fa_1 = pPCWAgg2d(1:11,dd_i==c);
    weights_fa_1 = weights_fa_1 ./ repmat(sum(weights_fa_1,2),1,size(weights_fa_1,2));
    weights_fa_2 = pPCWAgg2d(12:23,dd_i==c);
    weights_fa_2 = weights_fa_2 ./ repmat(sum(weights_fa_2,2),1,size(weights_fa_2,2));
    weights_fa_3 = pPCWAgg2d(24:35,dd_i==c);
    weights_fa_3 = weights_fa_3 ./ repmat(sum(weights_fa_3,2),1,size(weights_fa_3,2));
    weights_fa_4 = pPCWAgg2d(36:47,dd_i==c);
    weights_fa_4 = weights_fa_4 ./ repmat(sum(weights_fa_4,2),1,size(weights_fa_4,2));
    weights_fa_5 = pPCWAgg2d(48:end,dd_i==c);
    weights_fa_5 = weights_fa_5 ./ repmat(sum(weights_fa_5,2),1,size(weights_fa_5,2));

    % Average and write into prediction matrix.
    % First average over consumer types, then average over months within a
    % year.
    % Use of nanmean instead of mean is safety measure in case estimation
    % gest stuck in weird region where some consumer types / choices are not
    % simulated. This didn't happen with our final model specifications.
    mm_predict_fa_contract(0*8+c,:) = nanmean(...
                                           permute(sum(permute(ms_mm_aux_1,[1,3,2]) .* repmat(weights_fa_1,[1,1,J]),2),[1,3,2]) ... 
                                                ,1);
    mm_predict_fa_contract(1*8+c,:) = nanmean(...
                                           permute(sum(permute(ms_mm_aux_2,[1,3,2]) .* repmat(weights_fa_2,[1,1,J]),2),[1,3,2]) ... 
                                                ,1);
    mm_predict_fa_contract(2*8+c,:) = nanmean(...
                                           permute(sum(permute(ms_mm_aux_3,[1,3,2]) .* repmat(weights_fa_3,[1,1,J]),2),[1,3,2]) ... 
                                                ,1);
    mm_predict_fa_contract(3*8+c,:) = nanmean(...
                                           permute(sum(permute(ms_mm_aux_4,[1,3,2]) .* repmat(weights_fa_4,[1,1,J]),2),[1,3,2]) ... 
                                                ,1);
    mm_predict_fa_contract(4*8+c,:) = nanmean(...
                                           permute(sum(permute(ms_mm_aux_5,[1,3,2]) .* repmat(weights_fa_5,[1,1,J]),2),[1,3,2]) ... 
                                                ,1);
end

% (b) For partially informed consumers with contract type observed.
si_la_mm = cat(2,repmat(time_index,[1,1,NS_cons]), pSharesNoPCW3d);
mm_predict_la_contract = zeros(5*8,J);

% Average separately for each consumer type.
% Average over months within each year.
for c=1:8
    % Extract and separate predicted market shares by year.
    ms_mm_aux_1 = si_la_mm(1:11,2:end,dd_i==c,:); 
    ms_mm_aux_2 = si_la_mm(12:23,2:end,dd_i==c,:); 
    ms_mm_aux_3 = si_la_mm(24:35,2:end,dd_i==c,:); 
    ms_mm_aux_4 = si_la_mm(36:47,2:end,dd_i==c,:); 
    ms_mm_aux_5 = si_la_mm(48:end,2:end,dd_i==c,:); 

    % These weights indicate how important is an individual consumer
    % within the group of PCW non users.
    weights_la_1 = (1.0-pPCWAgg2d(1:11,dd_i==c));
    weights_la_1 = weights_la_1 ./ repmat(sum(weights_la_1,2),1,size(weights_la_1,2));
    weights_la_2 = (1.0-pPCWAgg2d(12:23,dd_i==c));
    weights_la_2 = weights_la_2 ./ repmat(sum(weights_la_2,2),1,size(weights_la_2,2));
    weights_la_3 = (1.0-pPCWAgg2d(24:35,dd_i==c));
    weights_la_3 = weights_la_3 ./ repmat(sum(weights_la_3,2),1,size(weights_la_3,2));
    weights_la_4 = (1.0-pPCWAgg2d(36:47,dd_i==c));
    weights_la_4 = weights_la_4 ./ repmat(sum(weights_la_4,2),1,size(weights_la_4,2));
    weights_la_5 = (1.0-pPCWAgg2d(48:end,dd_i==c));
    weights_la_5 = weights_la_5 ./ repmat(sum(weights_la_5,2),1,size(weights_la_5,2));

% Average and write into prediction matrix.
% First average over consumer types, then average over months within a
% year.
    mm_predict_la_contract(0*8+c,:) = nanmean(...
                                           permute(sum(permute(ms_mm_aux_1,[1,3,2]) .* repmat(weights_la_1,[1,1,J]),2),[1,3,2]) ... 
                                                ,1);
    mm_predict_la_contract(1*8+c,:) = nanmean(...
                                           permute(sum(permute(ms_mm_aux_2,[1,3,2]) .* repmat(weights_la_2,[1,1,J]),2),[1,3,2]) ... 
                                                ,1);
    mm_predict_la_contract(2*8+c,:) = nanmean(...
                                           permute(sum(permute(ms_mm_aux_3,[1,3,2]) .* repmat(weights_la_3,[1,1,J]),2),[1,3,2]) ... 
                                                ,1);
    mm_predict_la_contract(3*8+c,:) = nanmean(...
                                           permute(sum(permute(ms_mm_aux_4,[1,3,2]) .* repmat(weights_la_4,[1,1,J]),2),[1,3,2]) ... 
                                                ,1);
    mm_predict_la_contract(4*8+c,:) =  nanmean(...
                                           permute(sum(permute(ms_mm_aux_5,[1,3,2]) .* repmat(weights_la_5,[1,1,J]),2),[1,3,2]) ... 
                                                ,1);

end
% (c) For fully informed consumers with contract type unobserved.
% Recall that some firms only offer one contract.
mm_predict_fa_firm = zeros(5*8,7);
mm_predict_fa_firm(:,1) = mm_predict_fa_contract(:,1)+mm_predict_fa_contract(:,2);
mm_predict_fa_firm(:,2) = mm_predict_fa_contract(:,3)+mm_predict_fa_contract(:,4);
mm_predict_fa_firm(:,3) = mm_predict_fa_contract(:,5); 
mm_predict_fa_firm(:,4) = mm_predict_fa_contract(:,6)+mm_predict_fa_contract(:,7);
mm_predict_fa_firm(:,5) = mm_predict_fa_contract(:,8)+mm_predict_fa_contract(:,9);
mm_predict_fa_firm(:,6) = mm_predict_fa_contract(:,10);
mm_predict_fa_firm(:,7) = mm_predict_fa_contract(:,11);
% (d) For partially informed consumers with contract type unobserved.
mm_predict_la_firm = zeros(5*8,7);
mm_predict_la_firm(:,1) = mm_predict_la_contract(:,1)+mm_predict_la_contract(:,2);
mm_predict_la_firm(:,2) = mm_predict_la_contract(:,3)+mm_predict_la_contract(:,4);
mm_predict_la_firm(:,3) = mm_predict_la_contract(:,5);
mm_predict_la_firm(:,4) = mm_predict_la_contract(:,6)+mm_predict_la_contract(:,7);
mm_predict_la_firm(:,5) = mm_predict_la_contract(:,8)+mm_predict_la_contract(:,9);
mm_predict_la_firm(:,6) = mm_predict_la_contract(:,10);
mm_predict_la_firm(:,7) = mm_predict_la_contract(:,11);

%% Match model predictions to observed survey data.
% (a) For fully informed consumers with known contract type.
% Construct choice matrix.
% First column indicates consumer type (1-40) and second column
% indicates contract choice (1-11).
choice_full_known = dummyvar(mm_full_known(:,2))';
predict_full_known = mm_predict_fa_contract(mm_full_known(:,1),:)';
% Check whether individual choice probability predictions sum to one.
check_predict_full_known = abs(mean(sum(predict_full_known,1))-1.0);
if check_predict_full_known>1E-8
    fprintf('\n\nWARNING: Predicted choice probabilities for full-contract do not sum up to one!\n\n');
end

% (b) For fully informed consumers with unknown contract type.
% Construct choice matrix.
choice_full_unknown = dummyvar(mm_full_unknown(:,2))';
% Insert rows of zeros for firms for which we know contract type for
% sure.
choice_full_unknown = [choice_full_unknown(1:2,:); zeros(size(choice_full_unknown,2),1)'; choice_full_unknown(3:4,:); zeros(size(choice_full_unknown,2),1)'; choice_full_unknown(5,:)];
predict_full_unknown = mm_predict_fa_firm(mm_full_unknown(:,1),:)';
% Check that predicted choice probabilities sum to one.
check_predict_full_unknown = abs(mean(sum(predict_full_unknown,1))-1.0);
if check_predict_full_unknown>1E-8
    fprintf('\n\nWARNING: Predicted choice probabilities for full-firm do not sum up to one!\n\n');
end

% (c) For partially informed consumers with known contract type.
% Construct choice matrix.
choice_partial_known = dummyvar(mm_partial_known(:,2))';
predict_partial_known = mm_predict_la_contract(mm_partial_known(:,1),:)';
% Check that predicted choice probabilities sum to one.   
check_predict_partial_known = abs(mean(sum(predict_partial_known,1))-1.0);
if check_predict_partial_known>1E-8
    fprintf('\n\nWARNING: Predicted choice probabilities for partial-contract do not sum up to one!\n\n');
end

% (d) For partially informed consumers with unknown contract type.
% Construct choice matrix.
choice_partial_unknown = dummyvar(mm_partial_unknown(:,2))';
choice_partial_unknown = [choice_partial_unknown(1:2,:); zeros(size(choice_partial_unknown,2),1)'; choice_partial_unknown(3:4,:); zeros(size(choice_partial_unknown,2),1)'; choice_partial_unknown(5,:)];
predict_partial_unknown = mm_predict_la_firm(mm_partial_unknown(:,1),:)';
% Check that predicted choice probabilities sum to one.
check_predict_partial_unknown = abs(mean(sum(predict_partial_unknown,1))-1.0);
if check_predict_partial_unknown>1E-8
    fprintf('\n\nWARNING: Predicted choice probabilities for partial-firm do not sum up to one!\n\n');
end


%% Micromoments for PCW usage on individual level.
% Matrix for predicted PCW usage probabilities.
% years in rows, consumer types in columns.
mm_predict_active = zeros(5,8);
% Average separately for each consumer type, then average over months within each year.
    for c=1:8
        % This takes into account that PCW usage question is asked slightly differently over the different years.
        pcw_mm_aux_1 = pPCWAgg2d(1:7,dd_i==c); % Feb 2012- Sep 2012
        pcw_mm_aux_2 = pPCWAgg2d(1:19,dd_i==c); % Oct 2012 - Sep 2013
        pcw_mm_aux_3 = pPCWAgg2d(8:31,dd_i==c); % Oct 2013 - Sep 2014
        pcw_mm_aux_4 = pPCWAgg2d(20:43,dd_i==c); % Oct 2014 - Sep 2015
        % For last year, we cannot make predictions for July 2016 - Sep 2016 because of lack of aggregate data. Therefore, assume that the missing 3 months are not systematically different from the previous year and average over July 2015 - June 2016
        pcw_mm_aux_5 = pPCWAgg2d(30:end,dd_i==c); % Oct 2015 (-> July 2015) - Sep 2016 (-> June 2016)
        % Average (first over consumers in columns, then sum over months in year in rows)
        % Probability of having used PCW in a year is complement of not
        % having used PCW in any month of the year.
        % All adjustments are made in the lines above, so not ad-hoc scaling of probabilities should be done here.
        mm_predict_active(1,c) = (1-prod(1-mean(pcw_mm_aux_1,2)));
        mm_predict_active(2,c) = (1-prod(1-mean(pcw_mm_aux_2,2)));
        mm_predict_active(3,c) = (1-prod(1-mean(pcw_mm_aux_3,2)));
        mm_predict_active(4,c) = (1-prod(1-mean(pcw_mm_aux_4,2)));
        mm_predict_active(5,c) = (1-prod(1-mean(pcw_mm_aux_5,2)));
    end
% Here years are in rows, types in columns.
mm_predict_active = mm_predict_active';
% Just for diagnostics: Roughly weigh seniors and non-seniors 80%-20%
mm_predict_active_agg = 0.8*mean(mm_predict_active(1:4,:),1) + 0.2*mean(mm_predict_active(5:8,:),1);
% Match model predictions to observed data.      
pcw_use_predict = mm_predict_active(pcw_use_micro(:,1));
% PCW usage prediction error. PCW usage is recorded in column 4, i.e. end-1.
pcw_error_micro = pcw_use_predict - pcw_use_micro(:,end-1);

% Interact PCW use prediction error with regressors and instruments.
G_obs_4 = repmat(pcw_error_micro,1,n_mc_4) .* Z_4;
% END MM BLOCK: PCW Usage Moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MICROMOMENTS: Relative switching propensity on individual level.
% Matrix for predicted relative switching propensities.
% years in rows, consumer types in columns.
mm_predict_sw_ratio = zeros(5,8);
% Average churn rate across all consumer types for each survey period.    
sw_ratio_mean_1 = mean(mean(pChurnAgg2d(1:8,:)));
sw_ratio_mean_2 = mean(mean(pChurnAgg2d(9:20,:)));
sw_ratio_mean_3 = mean(mean(pChurnAgg2d(21:32,:)));
sw_ratio_mean_4 = mean(mean(pChurnAgg2d(33:46,:)));
sw_ratio_mean_5 = mean(mean(pChurnAgg2d(42:end,:)));
% Compute predicted churn rate for each consumer type (and year) separately.
for c=1:8
    % New version: survey period as relevant counterpart.
    sw_ratio_mm_aux_1 = pChurnAgg2d(1:8,dd_i==c);
    sw_ratio_mm_aux_2 = pChurnAgg2d(9:20,dd_i==c); 
    sw_ratio_mm_aux_3 = pChurnAgg2d(21:32,dd_i==c); 
    sw_ratio_mm_aux_4 = pChurnAgg2d(33:46,dd_i==c); 
    sw_ratio_mm_aux_5 = pChurnAgg2d(42:end,dd_i==c); 

    % Average (first over consumers in columns, then sum over months in year in rows)
    mm_predict_sw_ratio(1,c) = mean(mean(sw_ratio_mm_aux_1,2)) ./ sw_ratio_mean_1;
    mm_predict_sw_ratio(2,c) = mean(mean(sw_ratio_mm_aux_2,2)) ./ sw_ratio_mean_2;
    mm_predict_sw_ratio(3,c) = mean(mean(sw_ratio_mm_aux_3,2)) ./ sw_ratio_mean_3;
    mm_predict_sw_ratio(4,c) = mean(mean(sw_ratio_mm_aux_4,2)) ./ sw_ratio_mean_4;
    mm_predict_sw_ratio(5,c) = mean(mean(sw_ratio_mm_aux_5,2)) ./ sw_ratio_mean_5;
    end
% Transpose to have years are in rows, types in columns.
mm_predict_sw_ratio = mm_predict_sw_ratio';

% Match model predictions to observed data.  
% Extract relative switching ratio from micromoment data.    
% Note: relative switching propensity is coded in pcw_use_micro for easier
% code handling, even though this entails slightly confusing notation.
% This only extracts the consumer type index from survey data.
sw_ratio_mm_index = pcw_use_micro(:,1);
sw_ratio_predict = mm_predict_sw_ratio(sw_ratio_mm_index);
% Extract actual survey data on relative switching propensity.
sw_ratio_mm_data = pcw_use_micro(:,end);
% Relative switching propensity prediction error.
sw_ratio_error = sw_ratio_predict - sw_ratio_mm_data;
% Interact switching ratio prediction error with instruments.
G_obs_5 = repmat(sw_ratio_error,1,n_mc_5) .* Z_5;
%% END MM BLOCK: RELATIVE SWITCHING RATIO

% Construct G-matrices for micromoments on observation level.
% Firm and contract choices.
G_obs_mm{1} = (choice_full_known-predict_full_known)';
G_obs_mm{2} = (choice_full_unknown-predict_full_unknown)';
G_obs_mm{3} = (choice_partial_known-predict_partial_known)';
G_obs_mm{4} = (choice_partial_unknown-predict_partial_unknown)';
% END CONTRACT AND FIRM CHOICE MICROMOMENTS.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%x 
G_mm_1 = G_obs_mm{1}(:,1:end-1);
G_mm_2 = G_obs_mm{2}(:,1:end-1);
G_mm_3 = G_obs_mm{3}(:,1:end-1);
G_mm_4 = G_obs_mm{4}(:,1:end-1);

G_micro_61 = repmat(G_mm_1,1,n_mc_61) .* kron(Z_61,ones(1,J-1));
G_micro_62 = repmat(G_mm_2,1,n_mc_62) .* kron(Z_62,ones(1,7-1));
G_micro_63 = repmat(G_mm_3,1,n_mc_63) .* kron(Z_63,ones(1,J-1));
G_micro_64 = repmat(G_mm_4,1,n_mc_64) .* kron(Z_64,ones(1,7-1));

% Combine four different survey respondent groups into one array.
G_obs_6 = {G_micro_61,G_micro_62,G_micro_63,G_micro_64};

%% For easier diagnostics: Combine all moment matrices in one array.
G_obs_all = {G_obs_1,G_obs_2,G_obs_3,G_obs_4,G_obs_5,G_micro_61,G_micro_62,G_micro_63,G_micro_64};
for idx=1:length(G_obs_all)
    G_test = G_obs_all{idx};
    try 
        G_rank_test = rank(G_test)-size(G_test,2);
        if G_rank_test~=0
            fprintf('WARNING: Rank deficiency in moment set %d: %d\n',idx,G_rank_test);
        end
    catch
        warning('Moments contain NaNs. Assigning objective function value of 1E+12.');    
        % Load mean values from logit model.
        load(project_paths('OUT_ANALYSIS',['expmvalold_logit.mat']));
        save(project_paths('OUT_ANALYSIS',['expmvalold_',file_suffix,'.mat']),'expmvalold');
        % If expmvalold contain NaNs, exit function prematurely and reset delta values.
        model_myopic = 1E12;
        return;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute moment vector based on averaging over G matrices.
g_bar = [mean(G_obs_1,1), ... % BLP moments
         mean(G_obs_2,1), ... % aggregate churn rate moments 
         mean(G_obs_3,1), ... % aggregate PCW usage error moments
         mean(G_obs_4,1), ... % individual PCW usage error moments
         mean(G_obs_5,1), ... % individual switching heterogeneity moments
         mean(G_micro_61,1), ... % micro-level contract/firm choice error moments, group 1
         mean(G_micro_62,1), ... % micro-level contract/firm choice error moments, group 2
         mean(G_micro_63,1), ... % micro-level contract/firm choice error moments, group 3
         mean(G_micro_64,1), ... % micro-level contract/firm choice error moments, group 4
         ]; 

% Compute objective function value.
model_myopic = (g_bar * WMatrix * g_bar');
% Safety check for parameter guesses that result in NaN values for GMM function.
if isnan(model_myopic)==1
    model_myopic = 1E12;
    fprintf('\nObjective function not defined! Resetting Values!\n');
end
time = toc;
end