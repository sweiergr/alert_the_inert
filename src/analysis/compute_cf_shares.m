%{
    Compute counterfactual market shares when switching costs are
    contract-specific, i.e., work within a firm as well.
        % Required globals:
            - T, J, NS_cons, model_type
        % List of inputs:
            - optimal parameter vector: theta_opt_cf
            - preference shocks (for selecting averaging over different consumer types)
            - v_inc_cf_exp: matrix of exponentiated utilities for both
            - consumer-specific price coefficients.
            contracts for each consumer type.
        % List of outputs:
            - (Type-specific) Market shares for each product
            - Welfare matrix: averaged over all, average only over non-seniors,
            averaged only over seniors
%}
function [s_i, CS_inc_only_J]   = compute_cf_shares(theta_opt_cf,preference_shocks,v_inc_cf_exp, price_coeff_obs)
        
        global T NS_cons model_type
        % Define number of available products.
        J = size(v_inc_cf_exp,3);
        % Define vector of consumer-specific switching costs (depends on model
        % specification).
        if model_type==1 
            sc_i = repmat(theta_opt_cf(8),1,NS_cons);
        elseif model_type==4
            sc_i = (1-preference_shocks(1,:)) * theta_opt_cf(9) + preference_shocks(1,:) .* theta_opt_cf(10);
        elseif model_type==5
            sc_i = zeros(1,NS_cons);
        elseif model_type==6
            sc_i = repmat(theta_opt_cf(5),1,NS_cons);
        end
        % Prepare SC transformation to compute CCPs below.
        sc_i_exp = exp([zeros(NS_cons,1), sc_i', sc_i', zeros(NS_cons,1)]);
        % Container for consumer-specific market shares.
        s_ijt_cf_inc_only = zeros(J*T,NS_cons);
        % Set intial conditions for regulated market counterfactual (everybody
        % subsrcibe to default contract).
        s_ijt_0_reg = [ones(1,NS_cons);zeros(J-1,NS_cons)];
        s_lag = s_ijt_0_reg';
        % Init container for current period market share predictions.
        s_it = zeros(NS_cons,J);
        s_i = zeros(T,J,NS_cons);
        
        % Extract valuations (exponentiated) for both contracts.
        v_inc_cf_1 = squeeze(v_inc_cf_exp(:,:,1));
        v_inc_cf_1_green = squeeze(v_inc_cf_exp(:,:,2));
        
        CCP_t = zeros(J,J,NS_cons);
        % 3-d arrangement of switching cost terms.
        sc_3d = permute(repmat(sc_i',[,1,J,J]),[2,3,1]) .* repmat(~eye(J),[1,1,NS_cons]);
        sc_exp_3d = exp(sc_3d);
        % Loop over months.
        for t=1:T
            % Compute CCPs for current period and each consumer type.
            exp_v_t = repmat(squeeze(v_inc_cf_exp(t,:,:)),1,2);
            
            exp_v_it = permute(repmat(v_inc_cf_exp(1,:,:),[J,1,1]),[1,3,2]);
            
            CCP_t_num =  exp_v_it ./ sc_exp_3d;
            
            CCP_t_den = repmat(sum(CCP_t_num,2),[1,J,1]);
            
            
            CCP_t = CCP_t_num ./ CCP_t_den;
            for i =1:NS_cons
                s_it(i,:) = CCP_t(:,:,i)' * s_lag(i,:)';
            end
            
            % Reset s_lag for next period as current period market share
            % prediction.
            s_lag = s_it;
            % Write market share predictions for period t to total container.
            s_i(t,:,:) = s_it';
        end
        
        % Compute aggregate market share evolution.
        s_cf_inc_only = squeeze(mean(s_i,3));

        % Given predicted market shares, compute welfare. 
        % Now welfare depends on your previous contract.
        % If you are with contract 1 at beginning of period, welfare is given by:
        % Previous i-specific market shares.
        s_ijt_0 = reshape(s_ijt_0_reg,[1,J,500]);
        s_i_lag = cat(1,s_ijt_0,s_i(1:end-1,:,:));
        
        % Compute consumer surplus at each contract state and time period.
        CS_aux = zeros(T,J,NS_cons);
        CS_i = zeros(T,NS_cons);
        for t=1:T
            for j=1:J
                exp_v_t = squeeze(v_inc_cf_exp(t,:,:) ./ repmat(exp(sc_i),[1,1,J]));
                exp_v_t(:,j) = exp_v_t(:,j) .* exp(sc_i)';
                CS_aux(t,j,:) = log(sum(exp_v_t,2));
            end
            
            for i=1:NS_cons
                CS_i(t,i) = -(CS_aux(t,:,i) * s_i_lag(t,:,i)') ./ price_coeff_obs(i,1);
            end
        end

        % Finally, average over consumer types.
        CS_inc_only_J_all = mean(CS_i,2);
        CS_inc_only_J_young = mean(CS_i(:,preference_shocks(1,:)==0),2);
        CS_inc_only_J_old = mean(CS_i(:,preference_shocks(1,:)==1),2);
        
        % Combine welfare for different cases and demographic groups.
        CS_inc_only_J = [CS_inc_only_J_all, CS_inc_only_J_young, CS_inc_only_J_old];
end