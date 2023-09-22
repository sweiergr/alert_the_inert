%{ 

    This is the main file for the estimation of all model specifications reported in the paper. It loads and combines the various data sets and calls the file *model_myopic.m* to compute the value of GMM objective function. 

    This file generates the numbers for Table 2 in the main text, as well as figure 1.
    Formatting of the results is done in *src/final/format_est_results.py* and *src/final/plot_pcw_entrycost.py*.

%}
%% General preamble clearing memory and creating log-file.
clc
clear
clear global
%% Set global variables.
global  T J JF NS_cons NS_p_eps ...
        run_estimation SE_calc run_counterfactuals ...
        cf_subsidy_switch cf_perfect_info cf_reg_campaign cf_inc_only_all ...
        g_bar_all theta_all model_type debug_model

% Define list of models to estimate.
% For the paper we estimate models 1,4,5,6,7.
% 1: both psi and kappa are homogeneous.
% 2: psi is homogeneous, kappa is heterogeneous.
% 3: psi is heterogeneous, kappa is homogeneous.
% 4: both psi and kappa are heterogeneous.
% 5: only search frictions, i.e., no switching costs
% 6: only switching costs (capturing everything in reduced form), i.e., no search costs from using PCW.
% 7: baseline model without Hausman IV.
% 8: baseline model but using contraction mapping code with larger logit smoother.
model_list = [1,4,5,6,7,8];
% Define whether to two efficient two-step GMM or just first-step.
% rows indicate model to esitmate.
% first column indicates first step GMM, second column efficient two-step GMM.
gmm_type = [1,1; ... for baseline model 1
            1,0; ... for bigger model 2 (currently not used)
            1,0; ... for bigger model 3 (currently not used)
            1,1; ... for bigger model with both frictions heterogeneous across age groups
            ... in order to save some computation time, we only estimate first stage for the restricted models.
            ... second-stage estimates do not differ much and are available on request.
            1,0; ... for model without switching cost (strawman number 1)
            1,0; ... for model without PCW search costs (strawman number 2)
            1,1; ... for model without Hausman IV, otherwise identical to model 1.
            1,0; ... for baseline model (model 1) with larger logit smoother value (only for referee response).
            ];
% Set this to one if you want to see additional model diagnostics.
debug_model=0;
% Construct indicators for which system is used to estimate the model.
grid = 0; % this is only needed on machines where it's impossible to create the project_paths files.
% Extract the OS information
os_ind = computer;

if grid==1 % no waf: This is currently not needed anywhere since WY HPC runs also with project-path files.
    mex('contrMap_rcnl.c');
    code_path    = [pwd,'/'];
    results_path = [code_path(1:end-14),'/','bld/out/analysis/'];
    mkdir(results_path)
    logs_path    = [code_path(1:end-14),'/','bld/out/analysis/log/'];
    mkdir(logs_path)
    data_path    = [code_path(1:end-14),'/','bld/out/data/'];
else % i.e. everythign for which the waf-templates can be used.
    % This includes Wylie HPC after project_paths files have been build there manually.
    % Compile the C Code on the HPC Windows machines.
    % The mex compiling options may have to be adjusted for your machine.
    if os_ind(1:3)=='PCW' % Windows machines = Wylie HPC
        if isfolder('C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\')
            setenv('IFORT_COMPILER19','C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\');
            setenv('ICPP_COMPILER19','C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\');
        end
        if isfile('imsl_eval.dat')
            setenv('IMSL_LIC_FILE','imsl_eval.dat');
        end
       % This sets Windows computers up to use MS Visual C++ 2019 on WY HPC. 
       %mex -setup:C:\Users\sweiergr\AppData\Roaming\MathWorks\MATLAB\R2021a\mex_C_win64.xml C
       % This sets Windows computers up to use Intel Parallel Studio XE
       % 2019 to compile C code on WY HPC.
       %mex -setup:'C:\Program Files\MATLAB\R2021a\bin\win64\mexopts\intel_c_19_vs2019.xml' C
       % This compiles the C code with additional optimization options.
       mex('-R2018a','-I.\','COMPFLAGS="$COMPFLAGS"','OPTIMFLAGS="$OPTIMFLAGS /Ox"','LINKFLAGS="$LINKFLAGS"','LINKOPTIMFLAGS="$LINKOPTIMFLAGS /Ox"',project_paths('IN_ANALYSIS','contrMap_ES.c'))
       mex('-R2018a','-I.\','COMPFLAGS="$COMPFLAGS"','OPTIMFLAGS="$OPTIMFLAGS /Ox"','LINKFLAGS="$LINKFLAGS"','LINKOPTIMFLAGS="$LINKOPTIMFLAGS /Ox"',project_paths('IN_ANALYSIS','contrMap_ES_alt_LS.c'))
    elseif os_ind(1:3)=='GLN' % for Linux machines.
        mex('-R2018a',project_paths('IN_ANALYSIS','contrMap_ES.c'));
        mex('-R2018a',project_paths('IN_ANALYSIS','contrMap_ES_alt_LS.c'));
    elseif os_ind(1:3)=='MAC' % for Mac OS
        mex('-R2018a',project_paths('IN_ANALYSIS','contrMap_ES.c'));
        mex('-R2018a',project_paths('IN_ANALYSIS','contrMap_ES_alt_LS.c'));
    end
    % Safety measure: create new out-analysis folder.
    mkdir(project_paths('OUT_ANALYSIS','log'));
    code_path    = [pwd,'/'];
    results_path = project_paths('OUT_ANALYSIS');
    logs_path    = project_paths('OUT_ANALYSIS');
    data_path    = project_paths('OUT_DATA');
end

%% Set global parameters that determine what kind of estimation to run.
% Some basic diagnostics and statistics.
date = datestr(now);
print_details=1;
format 'short'

% Loop over models.
for model_idx=1:length(model_list)
    % Set which model type to run.
    model_type = model_list(model_idx);
    % Define model type file suffix.
    file_suffix = ['mod',num2str(model_type)];
    % Run GMM estimation with 2SLS weighting matrix block-diagonal.
    firststep_GMM = gmm_type(model_type,1);
    % Run GMM second stage using estimate of efficient weighting matrix based
    % on first step GMM results.
    efficient_GMM = gmm_type(model_type,2);
    % Set whether SE are computed for a variety of perturbations and
    % perturbation types.
    se_method = 'simple';
    % In order to make simulations exactly replicable, set RNG seed for simulated shocks.
    rng('default');
    % Globals to specify what to do.
    % Set this indicator to 1 if you want to do the estimation on your own
    % computer. If it is set to 0, it loads the estimated nonlinear parameters
    % directly from file.
    run_estimation = 1; % do not change this if you want to estimate the model yourself.
    SE_calc = 0; % switch to turn on parts in model code that are needed for SE computation.
    run_counterfactuals = 0; % switch to turn on parts in model code to run the counterfactuals.
    cf_subsidy_switch = 0; % indicator for switching subsidy counterfactual.
    cf_perfect_info = 0; % indicator for perfect information counterfactual.
    cf_reg_campaign = 0; % indicator for regulator campaign counterfactual.
    cf_inc_only_all = 0; % indicator for incumbent-all contracts counterfactual.

% Specify details fo which price IVs to use.
% Choose lags of weighted wholesale prices to use (in months).
whprice_lags = [3];
% Choose whether to use Hausman IV or not (0: no, 1: yes).
hausman_iv = 1;
if model_type==7
    hausman_iv = 0;
end

fprintf('Running model specification %d ...\n',model_type);
fprintf('Running first and second stage GMM? %d - %d \n', firststep_GMM, efficient_GMM);
%% Set path to diary log.
if grid==1 % currently not needed anywhere.
    input_path = 'input_data/';
    output_path = 'output_estimation/';
    output_path_log = output_path;
else
    % Add necessary paths to source libraries.
    input_path = project_paths('OUT_DATA');
    output_path = project_paths('OUT_TABLES');
    output_path_log = project_paths('OUT_ANALYSIS');
end
% Write log file to diary.
% Delete existing file if it already exists (to avoid appending).
if(exist([output_path_log,'estimation_log%02d']))
    delete([output_path_log,'estimation_log%02d']);
end
% Write new log-file.
diary([output_path_log,'estimation_log%02d']);
fprintf('This log contains background information on the estimation of discrete choice model with switching costs and endogenous PCW usage.\n')
fprintf(strcat('Model type estimated: %d.\n',model_type))
fprintf(strcat('Model was estimated on: \n',date))
fprintf('\nEstimation was computed on IU grid: %i \n', grid)

%% Load relevant data sets.
% Master data with market share data, prices, and churn rates.
fid = fopen([input_path,'new_master_data.csv'],'r');
% Manually adjust number of columns of data matrix: Currently 51.
C = textscan(fid, repmat('%s',1,29), 'delimiter',',', 'CollectOutput',true);
C = C{1};
fclose(fid);
% Transform to numerical data plus legend.
data = C(2:end,:);
data = cellfun(@str2double, data);

%% Define global model parameters and reshape elementary data components.
% Number of time periods/months = markets.
T = length(unique(data(:,1)));
% Number of products = contracts (constant across markets)
J = size(data,1) ./ T;
% Number of firms is hard-coded for now.
JF = 7;
% Number of simulated consumers.
NS_cons = 500;
% Number of simulation draws for price and epsilon shocks (simulation of
% PCW usage decision per consumer-market bin)
NS_p_eps = 30;
% Create legends for data matrices.
legend_main = C(1,:);
% Extract initial market shares
ms_initial = data(1:J,5);
% Reshape initial market shares.
s_init =  ms_initial';
% Extract market shares in wide format.
ms_aux = reshape(data(:,5),J,T)';
% Reshape observed market shares into column vector with outside good in
% there explicitly (needed for parts of C code).
market_shares = reshape(ms_aux',T*J,1);

% Initiate mean utility vector for first contraction mapping (at logit values).
delta = 1.0 * log( ms_aux(:,1:end-1) ./ repmat(ms_aux(:,end),1,J-1));
% Use this vector-representation of mean utilities
% Save exponentiated and add zero for outside good explicitly.
expmvalold = exp(reshape([delta, zeros(T,1)]',J*T,1));
% Delete mean utilities for first month.
expmvalold(1:J) = [];
% Save guess for exponentiated mean utilities to be loaded in model file.
save(project_paths('OUT_ANALYSIS',['expmvalold_logit.mat']),'expmvalold');
save(project_paths('OUT_ANALYSIS',['expmvalold_',file_suffix,'.mat']),'expmvalold');

%% Drop first period of data.
data(1:11,:) = [];
market_shares(1:J) = [];
T = T-1;
% Construct Nevo-style index vectors.
% Indicate to which market each observation belongs.
cdid = kron(linspace(1,T,T)',ones(J,1));
% Indicate index of last observation of each market.
cdindex = linspace(J,T*J,T)';
% Extract price beliefs in 100-EUR
price_beliefs = reshape(data(:,end),J,T)' ./100;
% Reshape churn rate data.
churn_rates = data(:,14);
new_month = zeros(length(data(:,1)),1);
new_month(:,1) = 1;
for t=2:length(data(:,1))
    new_month(t,1) = (data(t,1)~=data(t-1,1));
end
% Collapse to monthly aggregate churn rates.
avg_churn_rates = churn_rates(logical(new_month));
% END OF PREPARING MACRO DATA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load demographics data.
fid = fopen([input_path,'demo_data.csv'],'r');
% Manually adjust number of columns of data matrix: Currently 7.
CC = textscan(fid, repmat('%s',1,7), 'delimiter',',', 'CollectOutput',true);
CC = CC{1};
fclose(fid);
% Transform to numerical data plus legend.
legend_demo = CC(1,:);
data_demo = CC(2:end,:);
data_demo = cellfun(@str2double, data_demo);
% Drop first month of demographics data.
data_demo(1,:) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load advertising data.
fid = fopen([input_path,'advertisingdata.csv'],'r');
% Manually adjust number of columns of data matrix: Currently 5.
CC = textscan(fid, repmat('%s',1,5), 'delimiter',',', 'CollectOutput',true);
CC = CC{1};
fclose(fid);
% Transform to numerical data plus legend.
legend_adv = CC(1,:);
data_adv = CC(2:end,:);
data_adv = cellfun(@str2double, data_adv);
% Drop first month of advertising data (on firm level, not contract level).
data_adv(1:7,:) = [];
% Scale advertising data into millions.
data_adv(:,3) = data_adv(:,3) ./ 1000000;
data_adv(:,4) = data_adv(:,4) ./ 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load micromoments data from survey data.
fid = fopen([input_path,'micromoments.csv'],'r');
% Manually adjust number of columns of data matrix: Currently 8.
CC = textscan(fid, repmat('%s',1,8), 'delimiter',',', 'CollectOutput',true);
CC = CC{1};
fclose(fid);
% Transform to numerical data plus legend.
legend_mm = CC(1,:);
mm_data_raw = CC(2:end,:);
mm_data_raw = cellfun(@str2double, mm_data_raw);
% Count number of survey participants in each year.
micro_yu = unique(mm_data_raw(:,1));
survey_ny  = histc(mm_data_raw(:,1),micro_yu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data on initial conditions.
fid = fopen([input_path,'initialconditions.csv'],'r');
% Manually adjust number of columns of data matrix: Currently 10.
C = textscan(fid, repmat('%s',1,10), 'delimiter',',', 'CollectOutput',true);
C = C{1};
fclose(fid);
% Transform to numerical data plus legend.
legend_ic = C(1,:);
ic_data = C(2:end,:);
ic_data = cellfun(@str2double, ic_data(:,3:end));
% Labeling of columns:
% 1-4: nonseniors - 4 income groups as indicated by survey
% 5-8: seniors - 4 income groups as indicated by survey.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data on lagged wholesale prices.
fid = fopen([input_path,'ws_future_int_lags.csv'],'r');
% Manually adjust number of columns of data matrix: Currently 10.
C = textscan(fid, repmat('%s',1,10), 'delimiter',',', 'CollectOutput',true);
C = C{1};
fclose(fid);
% Transform to numerical data plus legend.
% Careful here: For some reason the cell array has an additional empty
% column.
wh_iv_lag_data = C(12:end,4:end);
wh_iv_lag_data = cellfun(@str2double, wh_iv_lag_data(:,1:end-1)) ./ 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load internet penetration and regulator campaign data.
fid = fopen([input_path,'internet_adv.csv'],'r');
% Manually adjust number of columns of data matrix: Currently 51.
CC = textscan(fid, repmat('%s',1,7), 'delimiter',',', 'CollectOutput',true);
CC = CC{1};
fclose(fid);
% Transform to numerical data plus legend.
legend_intadv = CC(1,:);
internet = CC(2:end,2:3);
internet = cellfun(@str2double, internet);
reg_active = CC(2:end,4);
reg_active = cellfun(@str2double, reg_active);
adv_agg = CC(2:end,5:6);
adv_agg = cellfun(@str2double, adv_agg)./10000000;
vtest_macro = CC(2:end,end);
vtest_macro = cellfun(@str2double, vtest_macro);
% Construct PCW use macro data.
% Here, we can expand set of moments/IVs used for PCW moments or specification changes.
% Specify how long effect of regulator campaign is.
% Baseline: It starts in month 9 (Oct 12 because campaign was started at end of Sep 12) and lasts for 6 months: elements 9:14 are one, everythign else is zero.
reg_active(1:8) = 0;
reg_active(9:14) = 1;
reg_active(15:end) = 0;
pcw_use_macro = [vtest_macro(2:end,:), ones(length(reg_active)-1,1), internet(2:end,:), reg_active(2:end,:), adv_agg(2:end,:)];
% END OF LOADING OF DATA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct and reshape regressor matrix.
% Currently this is just price and fixed effect for incumbent.
% 6: dummy for incumbent operator.
% 7: monthly subscription price (in 100-EUR)
% 8: dummy for green contract.
X_reg = data(:,6:8);
% Rescale price data.
X_reg(:,2) = X_reg(:,2) ./ 100;

% Rearrange regressors so that one can flexibly add firm fixed effects.
% Alternative specification that adds constant for inside goods to regressor matrix.
% instead of fixed effects.
% X_reg = [X_reg, repmat([ones(J-1,1);0],T,1)];
% Add firm fixed effects (instead of inside goods constant) to regressor matrix.
firm_fe = dummyvar(repmat([1;1;2;2;3;4;4;5;5;6;7],T,1));
contract_fe = dummyvar(repmat(linspace(1,J,J)',T,1));
% Adding time trend for incumbent.
time_trend = kron(linspace(1,T,T)'./10,ones(J,1));
time_trend_sq = time_trend.^2./10;
% Interact time trends with incumbent.
time_trend_inc = repmat(X_reg(:,1),1,2) .* [time_trend, time_trend_sq];
% For most recent version we estimate the model with firm fixed effects.
X_reg = [X_reg, firm_fe(:,2:end-1)];
% New order of regressors:
% 1. incumbent dummy: ECS-Electrabel
% 2. price -> not a linear parameter, so theta_linear is one element less.
% 3. green coefficient
% 4.-end firm FE or constant for inside goods.
% List of firm fixed effects.
% 4-1. EDF
% 5-2. Eneco (green only)
% 6-3. ENINuon
% 7-4. Essent
% 8-5. Lampiris (green only)

% Model with firm fixed effects: Initiate instrument matrix.
% Column 1: Dummy for incument.
% Column 2: green coefficient.
% Rest: Constant or firm fixed effects.
Z_matrix_BLP = X_reg(:,[1,3:end]);
% Construct price and switching cost instruments.
% This now takes a flexible argument for the number and type of lags.
% Set this at beginning of file.
price_iv = wh_iv_lag_data(:,whprice_lags);
% Append Hausman IV if desired.
if hausman_iv==1
    price_iv = [price_iv, data(:,12) ./1000];
end
sc_iv = data(:,23) ./ 100;
Z_matrix_BLP = [Z_matrix_BLP,  price_iv, sc_iv];


%% Construct 3-dimensional Halton draws for consumer simulation.
% Number of dimensions for Halton draws: NS_cons x N_rc + dim(adv_shocks) =
% 3 + 7*T = 3 + 7*53.
% Set the parameters for the Halton object.
HDRAWS = haltonset(374,'Skip',1000);
% Draw the random numbers for each simulated consumer.
halton_draws = net(HDRAWS,NS_cons);
% Normally distributed preference shock for green electricity.
% Invert uniform Halton draw to standard normal draw.
green_shocks = norminv(halton_draws(:,1)');
% Infer demographics from distribution: we assume that demogaphics do
% not change over time. This is very much consistent with the data that we
% have. 
% For age: 18.9% in market are seniors.
% Simulate seniors from second column of Halton draws.
age_demo = (halton_draws(:,2)'<0.189);
% Simulate income from third column of Halton draws.
income_shocks = halton_draws(:,3)';
% Actual mean and SD of Flemish income distribution in 1000-EUR (not in logs!).
income_mean = 2.5;
income_var = 2.25;
% Compute implied parameters mu and sigma of relevant lognormal
% distribution for income.
inc_mu = log(income_mean ./ sqrt(1+income_var./income_mean.^2));
inc_sigma2 = log(1 + income_var ./ income_mean.^2);
inc_sigma = sqrt(inc_sigma2);
% Create relevant lognormal distribution.
income_sim = logninv(income_shocks,inc_mu,inc_sigma);
% Precautionary measure: we top-code incoem at a monthly income of 15,000 EURs.
income_sim = min(income_sim,15.0);
% Sanity check: plot simulated income distribution.
hist(income_sim,75);
title('Simulated income distribution (100-EUR per month)')
saveas(gcf,project_paths('OUT_FIGURES','income_simulation.pdf'));
fprintf('Mean of simulated income distribution is %4.0f EUR.\n',1000*mean(income_sim));
fprintf('Median of simulated income distribution is %4.0f EUR.\n',1000*median(income_sim));
fprintf('Minimum of simulated income distribution is %4.0f EUR.\n',1000*min(income_sim));
fprintf('Maximum of simulated income distribution is %4.0f EUR.\n',1000*max(income_sim));
% Compute deviations from mean income and scale.
% This is deviation of income in levels, scaled in 100-EUR
dev_income = (income_sim - income_mean);
%% Combine preference shocks.
preference_shocks = [age_demo;dev_income; green_shocks;income_sim];
% Put demographic characteristics of simulatd consumers into table to
% select consumers for printing of CCP matrix.
char_label = {'c_idx','senior','income_dev','green_shock','income_sim'};
% Based on this, we pick the following representative consumers for illustrating CCP matrices
% Non-seniors: low-income-mean income-high income: 1-129-375
% Seniors: low-income-mean income-high income: 201-279-100
preference_table = table(linspace(1,NS_cons,NS_cons)',preference_shocks(1,:)', ...
    preference_shocks(2,:)',preference_shocks(3,:)',preference_shocks(4,:)', ...
    'VariableNames', char_label);
% Generate awareness-through-advertising shocks from last 7*53 columns of
% Halton draws.
awareness_shocks = reshape(halton_draws(:,4:end)',[T,7,NS_cons]);

% After preference shocks have been generated, expand initial conditions
% matrix.
% Initial conditions market share distribution is available for 8 consumer
% types -> assign each simulated consumer to one of the 8 types.
income_group = discretize(income_sim,[0,1.5,2.5,3.8,999999]);
income_group(income_group>4) = 4;
consumer_type = preference_shocks(1,:) * 4 + income_group;
if length(unique(consumer_type))<8
    fprintf('!!Warning!! Not all consumer types are simulated!!\nNumber of distinct consumer types simulated: %d\n.',length(unique(consumer_type)));
end
% Pick relevant column of initial market shares.
s_init_expand = ic_data(:,consumer_type);
% Into a long vector of J x I elements: first order is simulated consumers, then
% products.
s_init_long = reshape(s_init_expand,NS_cons*J,1);
% END OF SIMULATION OF CONSUMERS.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct advertising regressor matrix.
% Linear and squared terms of advertising expenditure.
adv_cons = ones(T,7,NS_cons);
% Several options: we use total monthly advertising.
% 3: actual monthly advertising.
% 4: normalized monthly advertising.
% 5: normalized monthly advertising per non-customer.
adv_X = repmat(reshape(data_adv(:,3),[7,T])',[1,1,NS_cons]);
adv_X_2 = adv_X.^2;
% Awareness effects due to age and income: being informed about market
% liberalization differently.
adv_X_senior = reshape(kron(preference_shocks(1,:),[ones(T,1),zeros(T,6)]),[T,7,NS_cons]);
adv_X_rich = reshape(kron(preference_shocks(2,:),[ones(T,1),zeros(T,6)]),[T,7,NS_cons]);
% Interaction terms between age and income and advertising expenditure.
adv_X_senior_int = adv_X_senior .* adv_X;
adv_X_rich_int = adv_X_rich .* adv_X;
% Keep model very simple in terms of advertising.
adv_data = cat(4, adv_cons, adv_X, adv_X_2, adv_X_senior, adv_X_rich, adv_X_senior_int, adv_X_rich_int);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct survey responeses for oberserved part of micromoments.
% Assign consumer type index to survey observations.
% Year-age-income
mm_inc_groups = 1 .* (mm_data_raw(:,4)<=3) + 2 .* (mm_data_raw(:,4) > 3 & mm_data_raw(:,4)<=5) ...
                + 3 .* (mm_data_raw(:,4) ==6) + 4 .* (mm_data_raw(:,4)>=7);
%% Create several aux index vectors.
% Consumer type-year index.
mm_ct_index = (mm_data_raw(:,1)-2012).*8 + mm_data_raw(:,5) .* 4 + mm_inc_groups;
% Consumer type index.
mm_c_index = mm_data_raw(:,5) .* 4 + mm_inc_groups;
% Year dummies.
mm_yrd = dummyvar(mm_data_raw(:,1)-2011);
% Consumer type dummies.
mm_c_type_dummies = dummyvar(mm_c_index);
mm_inc_dummies = dummyvar(mm_inc_groups);
mm_age_dummies = mm_data_raw(:,5);
% Prepare yearly average firm-level advertising as additional instrument
% for micro-level contract choice prediction error.
adv_exp_firm = adv_data(:,:,1,2);
adv_exp_contract = [adv_exp_firm(:,1), adv_exp_firm(:,1), ...
                        adv_exp_firm(:,2), adv_exp_firm(:,2), ...
                        adv_exp_firm(:,3), ...
                        adv_exp_firm(:,4), adv_exp_firm(:,4), ...
                        adv_exp_firm(:,5), adv_exp_firm(:,5), ...
                        adv_exp_firm(:,6), ...
                        adv_exp_firm(:,7) ];
% Containers for yearly statistics.
adv_yr_firm = zeros(5,7);
% Loop over years to aggregate monthly advertising expenditures.
for yr=1:5
    if yr==1
        adv_yr_firm(yr,:) = mean(adv_exp_firm(1:11,:));
    elseif yr>1 & yr<5 
        adv_yr_firm(yr,:) = mean(adv_exp_firm((yr-2)*12+12:yr*12-1,:));
    elseif yr==5
        adv_yr_firm(yr,:) = mean(adv_exp_firm((yr-2)*12+12:end,:));
    end
end
adv_yr_firm = adv_yr_firm';
mm_adv_idx = (mm_data_raw(:,1)-2012)*7+mm_data_raw(:,2);
mm_adv = adv_yr_firm(mm_adv_idx);
% Add instrument dummies to data matrix.
% We can experiment with different types of dummy IVs, i.p., income group dummies,
% age group dummies, or full consumer type dummies.
mm_data = [mm_data_raw, mm_ct_index,  mm_adv, mm_inc_dummies, mm_age_dummies];

% The coding is a bit counterintuitive, but more consistent with the data structure and makes it easier to extract information within the model code.
% 1: Consumer type index
% 2: Firm choice
% 3: Contract choice (if available)
% 4: vtest used (binary indicator)
% 5: relative switching propensity (to identify heterogeneity in switching
% cost)
pcw_use_micro = [mm_ct_index mm_data(:,2) mm_data(:,3) mm_data(:,6) mm_data(:,8)];

% Extract observed micromoments for PCW usage from survey data.
mm_data_active = zeros(8,5);
mm_data_weights = zeros(8,5);
mm_data_swhet = zeros(8,5);
for yr=1:5
    for ctype=1:8
        idx_select = (yr-1)*8+ctype;
        % Extract all consumers that fit characteristics.
        pcw_use_aux = pcw_use_micro(pcw_use_micro(:,1)==idx_select,end-1);
        swhet_aux = pcw_use_micro(pcw_use_micro(:,1)==idx_select,end);
        % Data on PCW usage.
        mm_data_active(ctype,yr) = mean(pcw_use_aux);
        % Data on switching heterogeneity.
        mm_data_swhet(ctype,yr) = mean(swhet_aux);
        mm_data_weights(ctype,yr) = length(pcw_use_aux);
    end
end
% Nomalize consumer type weights in survey for each year.
mm_data_weights = mm_data_weights ./ sum(mm_data_weights,1);
% Sanity check on observed data.
mm_data_active_agg = sum(mm_data_active .* mm_data_weights,1);
% Sanity check: How often do young people switch releative to old people?
mm_data_swhet_age_comp = mm_data_swhet(1:4,:) ./ mm_data_swhet(5:8,:);
% Export data to csv files for table formatting in Python.
dlmwrite([output_path_log,['mm_data_active_1_',file_suffix,'.csv']],mm_data_active,'delimiter',',');
dlmwrite([output_path_log,['mm_data_weights_1_',file_suffix,'.csv']],mm_data_weights,'delimiter',',');
dlmwrite([output_path_log,['mm_data_swhet_1_',file_suffix,'.csv']],mm_data_swhet,'delimiter',',');



% PCW instruments to interact with PCW use prediction error to form moment
% conditions. 
% - Demographics: dummies for income groups, dummies for senior
% - Dummy for choosing green contract.
% (- Dummy for choosing incumbent, firm dummies)
% (- Regulator campaign dummies)
% (- Internet penetration rate)
pcw_inst_micro = [dummyvar(mm_inc_groups) mm_data(:,5) (mm_data(:,3)==3) ];
% END OF PCW MACRO AND MICRO MOMENTS.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct micromoments based on firm/contract choice.

% Split survey data into four groups according to information &
% contract-type-known status.
% In addition, add relative switching frequency as additional micromoment
% data.
% New arrangement is: (1) consumer indicator (2) contract choice (3)
% relative switching frequency

% A. Consumers with known contract type and fully informed.
% Column 6 indicates whether consumer used vtest.
% Column 3 contains information about contract type. 0 indicates mising
% information on contract type.
mm_full_known = mm_data(mm_data(:,6)==1 & mm_data(:,3)~=0,:);
% Construct instruments for first group of survey participants.
% The Z_6i matrices contain the following instruments:
% - Firm-level advertising averaged for the year (1 column).
% - Consumer type dummies (8 columns)
% - Year dummies (5 dummies), last one will typically be dropped to avoid
% multicollinearity.
Z_6_full_known = mm_full_known(:,10:end);

mm_full_known = [mm_full_known, zeros(length(mm_full_known),1)];
% Construct combination of firm and contract type choice in one
% variable: firm * 10 + contract type.
mm_full_known(:,end) = mm_full_known(:,2).*10 + mm_full_known(:,3);
mm_full_known(mm_full_known(:,end)==11,end) = 1;
mm_full_known(mm_full_known(:,end)==13,end) = 2;
mm_full_known(mm_full_known(:,end)==21,end) = 3;
mm_full_known(mm_full_known(:,end)==23,end) = 4;
mm_full_known(mm_full_known(:,end)==33,end) = 5;
mm_full_known(mm_full_known(:,end)==41,end) = 6;
mm_full_known(mm_full_known(:,end)==43,end) = 7;
mm_full_known(mm_full_known(:,end)==51,end) = 8;
mm_full_known(mm_full_known(:,end)==53,end) = 9;
mm_full_known(mm_full_known(:,end)==63,end) = 10;
mm_full_known(mm_full_known(:,end)==71,end) = 11;
% Delete unnecessary columns: year, contract choices, demographics and
% vtest usage gets deleted
mm_full_known = [mm_full_known(:,9),mm_full_known(:,end), mm_full_known(:,8)];

% B. Consumers with unknown contract type and fully informed.
mm_full_unknown = mm_data(mm_data(:,6)==1 & mm_data(:,3)==0,:);
% Construct instruments for second group of survey participants.
Z_6_full_unknown = mm_full_unknown(:,10:end);
% This set of micromoments only requires consumer type information and
% firm choices and the relative switching frequency.
mm_full_unknown = [mm_full_unknown(:,9), mm_full_unknown(:,2),mm_full_unknown(:,8)]; 

% C. Consumers with known contract type and partially informed.
mm_partial_known = mm_data(mm_data(:,6)==0 & mm_data(:,3)~=0,:);
% Construct instruments for third group of survey participants.
Z_6_partial_known = mm_partial_known(:,10:end);

mm_partial_known = [mm_partial_known, zeros(length(mm_partial_known),1)];
mm_partial_known(:,end) = mm_partial_known(:,2).*10 + mm_partial_known(:,3);
mm_partial_known(mm_partial_known(:,end)==11,end) = 1;
mm_partial_known(mm_partial_known(:,end)==13,end) = 2;
mm_partial_known(mm_partial_known(:,end)==21,end) = 3;
mm_partial_known(mm_partial_known(:,end)==23,end) = 4;
mm_partial_known(mm_partial_known(:,end)==33,end) = 5;
mm_partial_known(mm_partial_known(:,end)==41,end) = 6;
mm_partial_known(mm_partial_known(:,end)==43,end) = 7;
mm_partial_known(mm_partial_known(:,end)==51,end) = 8;
mm_partial_known(mm_partial_known(:,end)==53,end) = 9;
mm_partial_known(mm_partial_known(:,end)==63,end) = 10;
mm_partial_known(mm_partial_known(:,end)==71,end) = 11;
% Arrange final micromoment data on consumer index, contract choice and
% relative switching frequency.
mm_partial_known = [mm_partial_known(:,9), mm_partial_known(:,end), mm_partial_known(:,8)];

% D. Consumers with unknown contract type and partially informed.
mm_partial_unknown = mm_data(mm_data(:,6)==0 & mm_data(:,3)==0,:);
% Construct instruments for fourth group of survey participants.
Z_6_partial_unknown = mm_partial_unknown(:,10:end);
mm_partial_unknown = [mm_partial_unknown(:,9), mm_partial_unknown(:,2),mm_partial_unknown(:,8)]; 

% Extract observed micro-level from cell computed in main file for 4
% different consumer types.
mm_data_reshaped = cell(4,1);
mm_data_reshaped{1} = mm_full_known;
mm_data_reshaped{2} = mm_full_unknown;
mm_data_reshaped{3} = mm_partial_known;
mm_data_reshaped{4} = mm_partial_unknown;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construction of objects needed for PCW usage simulation in model code.
% Draw logit shocks to simulate expected benefits of (not) searching
% PCW. Here order doesn't matter since everything is idd.
% Draw for each market, product, consumer type and simulation draw.
% Exponentiate iid shocks to speed up computation.
draws_eps = exp(evrnd(0,1,[T*J*NS_cons*NS_p_eps,1]));


%% Computation of price belief distribution parameters.
% The price belief distribution parameters are computed in
% merge_new_measterdata.do and saved in a separate Stata file.
% We allow the price belief distribution to change over
% time; variance is allowed to differ every 9 months.
% These numbers come from the auxiliary analysis of the price distribution in the Stata code.
sd_price_belief =  [0.04125814; ... empirical standard deviation of price distribution in first 9 months (in 100 EUR).
                   0.0331435; ...
                   0.0304667; ...
                   0.02335758; ...
                   0.02335235; ...
                   0.02726935];
month_per_period = (T+1) ./ length(sd_price_belief);
% Scale vector up to be consistent with mean of price beliefs.
% In main specification, we allow this to vary over 9-months periods, but not over
% firms.
sd_price_belief = kron(sd_price_belief,ones(month_per_period,1));
sd_price_belief(1,:) = [];
sd_price_belief = repmat(sd_price_belief, 1, J);

% For drawing from Gumbel distribution, we need to feed into evrnd:
% location parameter: mu and scale parameter: beta>0 
% We estimate: mean and SD from data and auxiliary regression.
% Translate these into EV-parameters by:
% 1. var = beta^2 pi^2 / 6
% Scale parameter is assumed to be constant over products:
beta_pb = sd_price_belief * sqrt(6) / pi;
% 2. mean = mu + beta * Euler-gamma
% Mean parameters are allowed to be different across contracts.
% Mean varies over products and months. 
mu_pb = price_beliefs - beta_pb * double(eulergamma);

% Version used for C code.
draws_p_normalized_long = zeros(T*J*NS_cons*NS_p_eps,1);
% Index vector for outside good.
idx_og = repmat([zeros(NS_p_eps*(J-1),1);ones(NS_p_eps,1)],NS_cons,1);

for t=1:T
    % Draw price beliefs from Gumbel distribution.
    % Location vector for period t resized to match efficient C code handling.
    % Order: consumer i (outer) - product k (middle) - simulation draw l (first).
    location_long_t = repmat(kron(mu_pb(t,:)',ones(NS_p_eps,1)),NS_cons,1);
    % Recall that in our main specification we allow the price distribution
    % variance to differ over time, but not over firms.
    draws_p_t_long = evrnd(location_long_t,beta_pb(t,1));
    outside_good_draws = draws_p_t_long(logical(idx_og));
    draws_outside_long = reshape(repmat(reshape(outside_good_draws,NS_p_eps,NS_cons),J,1),J*NS_cons*NS_p_eps,1);
    % Update price draws to reflect price differences between inside
    % and outside good price.
    draws_p_t_long_normalized = draws_p_t_long - draws_outside_long;
    % Simulate belief about price difference between all inside good and
    % outside good.
    draws_p_normalized_long((t-1)*J*NS_cons*NS_p_eps+1:t*J*NS_cons*NS_p_eps) = draws_p_t_long_normalized';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
end

price_beliefs = plot(draws_p_normalized_long);
title('Simulated normalized price belief realizations over time');
saveas(price_beliefs,project_paths('OUT_FIGURES',['price_beliefs_1_',file_suffix,'.pdf']));

% Use "long version" of data.
% Subtract actual price differences (inside  - outside good)!
% Format price data: be careful that effective price is price
% difference between outside good and inside good price.
% Substract outside good price here.
p_og_index = repmat([zeros(J-1,1);1],T,1);
p_og = X_reg(logical(p_og_index),2);
% Update column containing subscription price.
X_reg(:,2) = X_reg(:,2) - kron(p_og(:,1),ones(J,1));

% Difference price belief minus actual price and exponentiated draws.
% In model code, we multiply this difference with price coefficient and use
% it to simulate utility deviation due to price uncertainty/belief.
p_actual = kron(X_reg(:,2),ones(NS_cons*NS_p_eps,1));
draws_p_normalized = exp(draws_p_normalized_long-p_actual);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct instrument matrices for various sets of moments.
%{
    Block 1: BLP moments
    Block 2: Aggregate churn rates
    Block 3: Aggregate PCW usage
    Block 4: Micro/survey level PCW usage
    Block 5: Micro/survey based switching heterogeneity
    Block 6: Micro/suvey based contract choice errors (also interacted with advertising)
%}

% Create some auxiliary elements.
% Year dummies, quarter dummies, month-of-the-year dummies for macro data.
agg_year = [ones(11,1);kron([2;3;4],ones(12,1));5.*ones(6,1)];
aux_quarter = kron([1;2;3;4],ones(3,1));
aux_month = linspace(1,12,12)';
agg_quarter = [aux_quarter(2:end,1);repmat(aux_quarter,3,1);aux_quarter(1:6,1)];
agg_month =  [aux_month(2:end,1);repmat(aux_month,3,1);aux_month(1:6,1)];
% Respective dummy matrices.
agg_yd = dummyvar(agg_year);
agg_qd = dummyvar(agg_quarter);
agg_md = dummyvar(agg_month);
% Year dummies for micro data.
micro_yd = dummyvar(mm_data(:,1)-2011);

% 1. For BLP-xi - key element: instruments for prices and switching cost
Z_1 = Z_matrix_BLP;
% 2. For aggregate churn rate prediction error.
Z_2 = [agg_yd, agg_qd(:,1:end-1)];
% 3. For aggregate PCW usage error.
pcw_inst_macro = pcw_use_macro(:,[3,5,7]);
Z_3 = [agg_yd, agg_qd(:,1:end-1), pcw_inst_macro];
% 4. For micro PCW usage error.
Z_4 = [dummyvar(mm_c_index), micro_yd(:,2:end)];
% 5. For micro switching heterogeneity error.
Z_5 = [repmat(mm_data(:,8)>0,1,size(Z_4,2)) .* [dummyvar(mm_c_index), micro_yd(:,2:end)], ...
       repmat(mm_data(:,8)==0,1,size(Z_4,2)) .* [dummyvar(mm_c_index), micro_yd(:,2:end)]];
% 6. For micro contract choice prediction error
% Drop last year as reference year to avoid multicollinearity.
% Because some consumer types are not observed in the year 2015, we drop 2
% years from the micro moments IV matrix.
Z_61 = Z_6_full_known(:,1:end-1);
Z_62 = Z_6_full_unknown(:,1:end-1);
Z_63 = Z_6_partial_known(:,1:end-1);
Z_64 = Z_6_partial_unknown(:,1:end-1);
% Combine all instrument matrices in one array to pass to model file.
Z_matrix = {Z_1,Z_2,Z_3,Z_4,Z_5,Z_61,Z_62,Z_63,Z_64};
% Check instrument matrices for rank deficiencies.
rank_test_Z = zeros(size(Z_matrix,2),1);
for i=1:size(Z_matrix,2)
    Z_aux = Z_matrix{i};
    rank_test_Z(i) = size(Z_aux,2) - rank(Z_aux);
    if rank_test_Z(i)>0
        fprintf('Rank deficiency in IV matrix %d: %d linearly dependent columns.\n',i,rank_test_Z(i));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct weighting matrices.
% Number of survey respondents.
n_survey1 = size(mm_full_known,1);
n_survey2 = size(mm_full_unknown,1);
n_survey3 = size(mm_partial_known,1);
n_survey4 = size(mm_partial_unknown,1);
n_survey = n_survey1 + n_survey2 + n_survey3 + n_survey4;

% 1. Block: For xi-moments (BLP macro moments)
WM_21 = (J-1)*T*inv(Z_1'*Z_1);
% 2. Block: for aggregate churn rate error.
WM_22 = T * inv(Z_2'*Z_2);
% 3. Block: for aggregate PCW usage error.
WM_23 = T * inv(Z_3'*Z_3);
% 4. Block: for individual level PCW usage error.
WM_24 = n_survey * inv(Z_4'*Z_4);
% 5. Block: for relative switching ratio on individual level.
WM_25 = n_survey * inv(Z_5'*Z_5);
% 6. Block for individual level contract choice for different consumer
% groups (4 blocks).
WM_261 = n_survey1 * inv(Z_61'*Z_61);
WM_262 = n_survey2 * inv(Z_62'*Z_62);
WM_263 = n_survey3 * inv(Z_63'*Z_63);
WM_264 = n_survey4 * inv(Z_64'*Z_64);
% Concatenate block-diagonally for different contract prediction errors.
WM_261Cell = repmat({WM_261}, 1, J-1);
WM_261J = blkdiag(WM_261Cell{:});
WM_262Cell = repmat({WM_262}, 1, 7-1);
WM_262J = blkdiag(WM_262Cell{:});
WM_263Cell = repmat({WM_263}, 1, J-1);
WM_263J = blkdiag(WM_263Cell{:});
WM_264Cell = repmat({WM_264}, 1, 7-1);
WM_264J = blkdiag(WM_264Cell{:});
WM_26 = blkdiag(WM_261J,WM_262J,WM_263J,WM_264J);

% Extract matrix of advertising expenditures to construct additional
% micromoments on yearly level to match with micromoment frequency.
% Construct index matrices to select correct advertising amount.
yr_id_full_known = zeros(size(mm_full_known,1),J);
yr_id_full_unknown = zeros(size(mm_full_unknown,1),7);
yr_id_partial_known = zeros(size(mm_partial_known,1),J);
yr_id_partial_unknown = zeros(size(mm_partial_unknown,1),7);
for j=1:J
    yr_id_full_known(:,j) = (j-1)*5+ceil(mm_full_known(:,1)./8);
    yr_id_partial_known(:,j) = (j-1)*5+ceil(mm_partial_known(:,1)./8);
end
for frm=1:7
    yr_id_full_unknown(:,frm) = (frm-1)*5+ceil(mm_full_unknown(:,1)./8);
    yr_id_partial_unknown(:,frm) = (frm-1)*5+ceil(mm_partial_unknown(:,1)./8);
end

adv_exp_firm = adv_data(:,:,1,2);
adv_exp_contract = [adv_exp_firm(:,1), adv_exp_firm(:,1), ...
                        adv_exp_firm(:,2), adv_exp_firm(:,2), ...
                        adv_exp_firm(:,3), ...
                        adv_exp_firm(:,4), adv_exp_firm(:,4), ...
                        adv_exp_firm(:,5), adv_exp_firm(:,5), ...
                        adv_exp_firm(:,6), ...
                        adv_exp_firm(:,7) ];
% Containers for yearly statistics.
adv_yr_firm = zeros(5,7);
adv_yr_contract = zeros(5,J);
% Loop over years to aggregate monthly advertising expenditures.
for yr=1:5
    if yr==1
        adv_yr_firm(yr,:) = mean(adv_exp_firm(1:11,:));
        adv_yr_contract(yr,:) = mean(adv_exp_contract(1:11,:));

    elseif yr>1 & yr<5 
        adv_yr_firm(yr,:) = mean(adv_exp_firm((yr-2)*12+12:yr*12-1,:));
        adv_yr_contract(yr,:) = mean(adv_exp_contract((yr-2)*12+12:yr*12-1,:));

    elseif yr==5
        adv_yr_firm(yr,:) = mean(adv_exp_firm((yr-2)*12+12:end,:));
        adv_yr_contract(yr,:) = mean(adv_exp_contract((yr-2)*12+12:end,:));
    end
end
adv_full_known = adv_yr_contract(yr_id_full_known);
adv_full_unknown = adv_yr_firm(yr_id_full_unknown);
adv_partial_known = adv_yr_contract(yr_id_partial_known);
adv_partial_unknown = adv_yr_firm(yr_id_partial_unknown);
% Combine 4 groups of advertising data into one object.
adv_data_mm = cell(4,1);
adv_data_mm{1} = adv_full_known;
adv_data_mm{2} = adv_full_unknown;
adv_data_mm{3} = adv_partial_known;
adv_data_mm{4} = adv_partial_unknown;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stack in block diagonal matrix.
WM_2sls = blkdiag(WM_21, WM_22, WM_23, WM_24, WM_25, WM_26);


%% Supply weighting matrix directly to objective function.
% 1 for identity matrix.
% 2 for 2SLS (block diagonal)
% 3 for efficient WM (requires first-stage based on 1 or 2 + second run
WMatrix = WM_2sls;
WM2sls_CN = cond(WMatrix);
fprintf('Condition number of first-step weighting matrix is %e.\n',WM2sls_CN);
fprintf('Is 2SLS weighting matrix rank deficient? %d linearly dependent columns!\n',rank(WMatrix)-size(WMatrix,2));

%% Diagnostics of initial conditions.
ms_wide = reshape(market_shares,J,T);
s_init_db_micro = reshape(s_init_long,J,NS_cons);
s_init_db_agg =  repmat(s_init',1,NS_cons);
% Compare implied aggregate initial conditions based on micro data.
s_init_agg_micro = mean(reshape(s_init_long,J,NS_cons),2);
s_init_agg_agg = s_init';
% Compute weights based on micro level initial conditions.
s_init_weights = s_init_db_micro ./ repmat(s_init_agg_micro,1,NS_cons);
s_init_weights(isnan(s_init_weights)==1) = 0;
s_init_weighted = s_init_db_agg .* s_init_weights;
s_init_weighted_norm = s_init_weighted ./ repmat(sum(s_init_weighted,1),J,1);
% Use this for now: initial conditions that take into account relative
% shares from survey but normalize such that they are very close to
% aggregate initial market shares.
s_init_long = reshape(s_init_weighted_norm,NS_cons*J,1);

%{ 
    Load starting values for first step GMM estimation from corresponding files in src/model_specs.
    Order of parameters for model 1 and 7:
    1-3: PCW usage costs
    4. mean price coefficient (nonlinear)
    5. standard deviation of green preference
    6. incumbent-age interaction
    7. price-income interaction
    8. switching cost
    9-10: advertising process parameters

    Order of parameters for model 4:
    1-4: PCW usage costs
    5. mean price coefficient (nonlinear)
    6. standard deviation of green preference
    7. incumbent-age interaction
    8. price-income interaction
    9-10. switching cost
    11-12: advertising process parameters

    Order for models 5 identical to one except that the restricted parameters are omitted.
%}

load(project_paths('IN_MODEL_SPECS',['theta0_',file_suffix]));

% Specify convergence criteria for contraction mapping routine and define
% anoynmous function.
% Here we just solve the model a few times to explore the properties of the contraction mapping.
tol_delta_c = 1E-6;
max_iter_c = 650;
gmm_function = @(theta) model_myopic(theta, cdid, cdindex, market_shares, s_init_long, avg_churn_rates, mm_data_reshaped, pcw_use_micro, pcw_inst_micro, pcw_use_macro, adv_data, adv_data_mm, data_demo, X_reg, Z_matrix, draws_p_normalized, draws_eps,  preference_shocks,awareness_shocks,WMatrix, tol_delta_c,max_iter_c);
% Evaluate model at starting values.
gmm_obj_0 = gmm_function(theta_0);
% Load options for actual model estimation.
load(project_paths('IN_MODEL_SPECS',['c_options_',file_suffix,'1']));
gmm_function = @(theta) model_myopic(theta, cdid, cdindex, market_shares, s_init_long, avg_churn_rates, mm_data_reshaped, pcw_use_micro, pcw_inst_micro, pcw_use_macro, adv_data, adv_data_mm, data_demo, X_reg, Z_matrix, draws_p_normalized, draws_eps,  preference_shocks,awareness_shocks,WMatrix, tol_delta_c,max_iter_c);
[~,~,~,~,~,~,~, ...
pred_shares1d_1, pred_churnAgg1d_1, ...
pred_cSurplusPCW_1,pred_cSurplusNoPCW_1, ...
pred_PCWAgg1d_1,~,~, ...
pred_CCPAvg_1] = gmm_function(theta_0);
    
%% Evaluate aggregate vtest numbers (just for initial diagnostics).
churn_wide = reshape(churn_rates,J,T)';
churn_ts = churn_wide(:,1);
month_grid = linspace(1,T,T)';
vtest_check = [ month_grid, pcw_use_macro(:,1), churn_ts, pcw_use_macro(:,5)];
plot(month_grid(2:end),pcw_use_macro(2:end,1),month_grid(2:end),churn_ts(1:end-1));
legend('PCW Usage','Churn Rate');
title('Comparison of PCW usage rate and churn rate');

%% Initialize diagnostics matrices for evaluated parameter values.
% Write initial guesses for delta to matrix that collects all the mean
% utilities for each iteration.
theta_all = theta_0;
g_bar_all = zeros(1,size(WMatrix,2));

%% Call optimizer.
% Run first step GMM.
if firststep_GMM==1
    if run_estimation==1
        
    theta_opt = theta_0;
    theta_opt_ars = theta_opt;
    % (1) fminunc is gradient-based, therefore very efficient, but can also be
    % tricky if objective function doesn't have a well behaved gradient.
    %[theta_opt, gmm_value] = fminsearch(gmm_function, theta_0,fminsearch_options);
    % Run model file to update linear parameters.
    %gmm_function(theta_opt);
    
    % (2) ARS (accelerated random search) is also a line-search algorithm, but
    % more sophisticated than pure simplex search. We used this to get a rough idea of the potential for local minima 
    % across a broad search radius.
    % Set parameter bounds and options for ARS optimization routine.
    % We use this to roughly check a broad range of parameter values and
    % then use fminsearch to verify that the ARS minimum is indeed a
    % minimum of the objective function.
    % For debugging and exploration of parameter space: 
    LB = theta_0 - 3.75;
    UB = theta_0 + 3.75;
    % Just for diagnostics: check where optimizers starts.   
    theta_ars_start = 0.5 .* LB + 0.5 .* UB;
    % Set contraction parameter for ARS.
    c = 2;
    % Set convergence criterion for ARS.
    rho = 10^-6;
    %% Actual minimization starts is done here.
    %theta_opt = ars(gmm_function,LB,UB, c,rho,1000,5);
    % Transpose because ars.m outputs row vector.
    % theta_opt = theta_opt';
    % theta_opt_ars = theta_opt;
    % param_check_ars_1 = [LB, theta_0, theta_opt_ars, UB];
    
    % (2) fminsearch is line-search algorithm, that is relatively slow, but more
    % robust than gradient-based optimization, so probably the better choice.
    % Load fminsearch options.
    load(project_paths('IN_MODEL_SPECS',['fminsearch_options_',file_suffix,'1']));
    format long
    % Start actual first-stage GMM estimation.
    [theta_opt, gmm_value, opt_flag_1, opt_out_1] = fminsearch(gmm_function, theta_opt,fminsearch_options);
    format short
    % Save nonlinear parameters after first-stage estimation.
    dlmwrite([output_path_log,['theta_opt_1_nl_',file_suffix,'.csv']],theta_opt,'delimiter',',');
    % Save linear parameters after first-stage estimation.
    load(project_paths('OUT_ANALYSIS',['thetalinear_',file_suffix,'.mat']),'theta_linear');
    theta_lin_opt_1 = theta_linear;
    save(project_paths('OUT_ANALYSIS',['theta_lin_opt_1_',file_suffix]),'theta_lin_opt_1');
    % Save mean utilities after first-stage estimation.
    load(project_paths('OUT_ANALYSIS',['expmvalold_',file_suffix,'.mat']),'expmvalold');
    expmval_opt_1 = expmvalold;
    save(project_paths('OUT_ANALYSIS',['expmvalold_opt_1_',file_suffix]),'expmval_opt_1');
        
    %% Some diagnostics in case of convergence issues.
    pshares_track = [market_shares];
    cSurplusPCW_track = [];
    cSurplusNoPCW_track = [];
    expmvalold_track = expmval_opt_1;
    theta_linear_track = theta_lin_opt_1;
    % Do another model evaluation at optimal parameter values to ensure convergence
    % and get predicted market shares, churn rates, and PCW usage.
    % You shouldn't see any changes here.
    tol_delta_c = 1E-4;
    max_iter_c = 750;
        gmm_function = @(theta) model_myopic(theta, cdid, cdindex, market_shares, s_init_long, avg_churn_rates, mm_data_reshaped, pcw_use_micro, pcw_inst_micro, pcw_use_macro, adv_data, adv_data_mm, data_demo, X_reg, Z_matrix, draws_p_normalized, draws_eps,  preference_shocks,awareness_shocks,WMatrix, tol_delta_c,max_iter_c);
        [~,~,~,~,~,~,~, ...
        pred_shares1d_1, pred_churnAgg1d_1, ...
        pred_cSurplusPCW_1,pred_cSurplusNoPCW_1, ...
        pred_PCWAgg1d_1,~,~, ...
        pred_CCPAvg_1] = gmm_function(theta_opt);
    
    % Check that linear parameters and mean utilities have not changed.
    load(project_paths('OUT_ANALYSIS',['thetalinear_',file_suffix,'.mat']),'theta_linear');
    load(project_paths('OUT_ANALYSIS',['expmvalold_',file_suffix,'.mat']),'expmvalold');
    conv_check_meanutil_1 = norm(expmval_opt_1-expmvalold);
    conv_check_thetalin_1 = norm(theta_lin_opt_1-theta_linear);
    
    % Track changes to various model statistics.
    pshares_track = [pshares_track, pred_shares1d_1];
    expmvalold_track = [expmvalold_track,expmvalold];
    theta_linear_track = [theta_linear_track, theta_linear];
    % Check changes across iterations.
    change_pshares = pshares_track(:,end) ./ pshares_track(:,end-1);
    change_expmval = expmvalold_track(:,end) ./ expmvalold_track(:,end-1);
    change_thetalinear = theta_linear_track(:,end) ./ theta_linear_track(:,end-1);
    
    % Do one final  model evaluation with a somewhat tighter convergence criterion at optimal parameter values to ensure convergence during estimation. You shouldn't see any changes here either.
    tol_delta_c = 1E-8;
    max_iter_c = 1000;
    gmm_function = @(theta) model_myopic(theta, cdid, cdindex, market_shares, s_init_long, avg_churn_rates, mm_data_reshaped, pcw_use_micro, pcw_inst_micro, pcw_use_macro, adv_data, adv_data_mm, data_demo, X_reg, Z_matrix, draws_p_normalized, draws_eps,  preference_shocks,awareness_shocks,WMatrix, tol_delta_c,max_iter_c);
    [~,~,~,~,~,~,~, ...
    pred_shares1d_1, pred_churnAgg1d_1, ...
    pred_cSurplusPCW_1,pred_cSurplusNoPCW_1, ...
    pred_PCWAgg1d_1,~,~, ...
    pred_CCPAvg_1,mm_pred_active_1,mm_pred_swhet_1] = gmm_function(theta_opt);
    
    % Export predicted micromoments for PCW usage and switching ratios to csv files for table formatting in Python.
    dlmwrite([output_path_log,['mm_pred_active_1_',file_suffix,'.csv']],mm_pred_active_1,'delimiter',',');
    dlmwrite([output_path_log,['mm_pred_swhet_1_',file_suffix,'.csv']],mm_pred_swhet_1,'delimiter',',');
    % Aggregate predictions to the year level.
    mm_pred_active_1_agg = sum(mm_data_weights .* mm_pred_active_1,1);
    % Compute prediction difference.
    mm_data_active_diff = mm_pred_active_1_agg -mm_data_active_agg;
    
    % Compute yearly share of pcw users.
    % Assuming each household only uses it once.
    pcw_usage_agg_yr_data = [sum(pcw_use_macro(1:7,1)), sum(pcw_use_macro(1:19,1)), ...
                        sum(pcw_use_macro(8:31,1)), sum(pcw_use_macro(20:43,1)), ...
                        sum(pcw_use_macro(30:end,1))];
    pcw_usage_agg_yr_pred = [sum(pred_PCWAgg1d_1(1:7,1)), sum(pred_PCWAgg1d_1(1:19,1)), ...
                        sum(pred_PCWAgg1d_1(8:31,1)), sum(pred_PCWAgg1d_1(20:43,1)), ...
                        sum(pred_PCWAgg1d_1(30:end,1))];
    % Compare macro and micro level PCW usage rates.
    % First two rows: observed macro PCW usage and micro-level PCW usage.
    % Second two rows: predicted versions of the same.
    compare_macro_micro_pcw = [pcw_usage_agg_yr_data; mm_data_active_agg; ... 
                               pcw_usage_agg_yr_pred; mm_pred_active_1_agg];
    
    compare_macro_micro_pcw_mtly = compare_macro_micro_pcw ./ [8,12,12,12,12];
    
    compare_macro_micro_diff = (compare_macro_micro_pcw(3:4,:) - compare_macro_micro_pcw(1:2,:))./ [9,12,12,12,12];
   
    % Check that linear parameters and mean utilities have not changed.
    load(project_paths('OUT_ANALYSIS',['thetalinear_',file_suffix,'.mat']),'theta_linear');
    load(project_paths('OUT_ANALYSIS',['expmvalold_',file_suffix,'.mat']),'expmvalold');
    conv_check_meanutil_1 = norm(expmval_opt_1-expmvalold);
    conv_check_thetalin_1 = norm(theta_lin_opt_1-theta_linear);
    % Track changes to various model statistics.
    pshares_track = [pshares_track, pred_shares1d_1];
%     cSurplusPCW_track = [cSurplusPCW_track, pred_cSurplusPCW_1];
%     cSurplusNoPCW_track = [cSurplusNoPCW_track, pred_cSurplusNoPCW_1];
    expmvalold_track = [expmvalold_track,expmvalold];
    theta_linear_track = [theta_linear_track, theta_linear];
    
    % Check changes across iterations.
    change_pshares = pshares_track(:,end) ./ pshares_track(:,end-1);
%     change_csurpluspcw = cSurplusPCW_track(:,end) ./ cSurplusPCW_track(:,end-1);
%     change_csurplusnopcw = cSurplusNoPCW_track(:,end) ./ cSurplusNoPCW_track(:,end-1);
    change_expmval = expmvalold_track(:,end) ./ expmvalold_track(:,end-1);
    change_thetalinear = theta_linear_track(:,end) ./ theta_linear_track(:,end-1);
    % Final save nonlinear parameters after first-stage estimation.
    dlmwrite([output_path_log,['theta_opt_1_nl_',file_suffix,'.csv']],theta_opt,'delimiter',',');
    % Save linear parameters after first-stage estimation.
    load(project_paths('OUT_ANALYSIS',['thetalinear_',file_suffix,'.mat']),'theta_linear');
    theta_lin_opt_1 = theta_linear;
    save(project_paths('OUT_ANALYSIS',['theta_lin_opt_1_',file_suffix]),'theta_lin_opt_1');
    % Save mean utilities after first-stage estimation.
    load(project_paths('OUT_ANALYSIS',['expmvalold_',file_suffix,'.mat']),'expmvalold');
    expmval_opt_1 = expmvalold;
    save(project_paths('OUT_ANALYSIS',['expmvalold_opt_1_',file_suffix]),'expmval_opt_1');
    
%% Extract coefficients.
% Mean price coefficient and interaction with income.
% These depend on model specification.
% For baseline model (1), model without Hausman IV (7), and baseline with larger logit smoother (8)
if model_type==1 | model_type==7 || model_type==8
    theta_kappa = theta_opt(1:3);
    theta_price = [theta_opt(4); theta_opt(7)];
    theta_incumbent = [theta_linear(1); theta_linear(1) + theta_opt(6)];
    theta_green = [theta_linear(2); abs(theta_opt(5))];
    theta_sc = theta_opt(8);
    theta_adv = theta_opt(9:end);
% Not used in final version.
elseif model_type==2
    theta_kappa = theta_opt(1:4);
    theta_price = [theta_opt(5); theta_opt(8)];
    theta_incumbent = [theta_linear(1); theta_linear(1) + theta_opt(7)];
    theta_green = [theta_linear(2); abs(theta_opt(6))];
    theta_sc = theta_opt(9);
    theta_adv = theta_opt(10:end);
% Not used in final version.
elseif model_type==3
    theta_kappa = theta_opt(1:3);
    theta_price = [theta_opt(4); theta_opt(7)];
    theta_incumbent = [theta_linear(1); theta_linear(1) + theta_opt(6)];
    theta_green = [theta_linear(2); abs(theta_opt(5))];
    theta_sc = theta_opt(8:9);
    theta_adv = theta_opt(10:end);
% Model extension with heterogeneous market frictions.
elseif model_type==4
    theta_kappa = theta_opt(1:4);
    theta_price = [theta_opt(5); theta_opt(8)];
    theta_incumbent = [theta_linear(1); theta_linear(1) + theta_opt(7)];
    theta_green = [theta_linear(2); abs(theta_opt(6))];
    theta_sc = theta_opt(9:10);
    theta_adv = theta_opt(11:end);
% Model without switching costs.
elseif model_type==5
    theta_kappa = theta_opt(1:3);
    theta_price = [theta_opt(4); theta_opt(7)];
    theta_incumbent = [theta_linear(1); theta_linear(1) + theta_opt(6)];
    theta_green = [theta_linear(2); abs(theta_opt(5))];
    theta_sc = 0;
    theta_adv = theta_opt(8:end);
% Model without PCW search costs.
elseif model_type==6
    theta_kappa = zeros(3,1);
    theta_price = [theta_opt(1); theta_opt(4)];
    theta_incumbent = [theta_linear(1); theta_linear(1) + theta_opt(3)];
    theta_green = [theta_linear(2); abs(theta_opt(2))];
    theta_sc = theta_opt(5);
    theta_adv = theta_opt(6:end);
end

%% Print results to screen.
    fprintf('Estimates for price coefficients for mean price and income-interaction coefficient: \n %6.4f \n %6.4f \n %6.4f \n', ...
            theta_price(1),theta_price(2))
    fprintf('Estimates for incumbent preference for young and senior: \n %6.4f \n %6.4f \n', ...
            theta_incumbent(1), theta_incumbent(2))
    fprintf('Estimates for green coefficients (mean and variance): \n !! CAREFUL WITH INTERPRETATION OF THIS IF PREFERENCE IS LOGNORMAL!!\n \n %6.4f \n %6.4f \n', ...
            theta_green(1), theta_green(2))
    if model_type==1 || model_type==2 || model_type==7 || model_type==8
        fprintf('Estimates for switching costs: \n %6.4f \n', ...
                theta_sc(1))
    elseif model_type==3 || model_type==4
        fprintf('Estimates for switching costs (young and old, respectively): \n %6.4f \n %6.4f \n', ...
                theta_sc(1), theta_sc(2))
    end
    fprintf('Estimates for advertising coefficients: \n %6.4f \n %6.4f \n', ...
             theta_adv(1), theta_adv(2))

% Calculate implied monetary magnitudes (for quantity model with log prices).
% Willingness to pay for incumbent.
inc_EUR = -100 .* theta_incumbent ./ repmat(theta_price(1),2,1);
% Switching costs in EUR.
sc_EUR = -100 * theta_sc ./ theta_price(1);
% Willingness to pay for green electricity (mean).
% Auxiliary computations for specification where deviation if modeled as
% lognormal and total preference is given by "mean preference" (=lower
% bound, theta_green(1)) + lognormal deviation with parameters 0 and theta_green(2).
% Compute mean and variance of consumer-specific upward preference deviation.
green_dev_mean = exp(0.5 * theta_green(2)^2);
green_dev_var = exp( 2* theta_green(2)^2 - 1);
% If green preference is specified as log-normal: Mean preference for green electricity.
green_EUR = -100 ./ theta_price(1) .* (theta_green(1) + green_dev_mean);
green_lower_EUR = -100 ./ theta_price(1) .* theta_green(1);
% Only if green preference is specified as normal.
% green_EUR = -100 .* theta_green(1) ./ theta_price(1);

% Compute distribution of green electricity coefficient.
% New version where green preference is modeled as log-normal distribution.
gc_grid = linspace(0,15,100)';
gc_dist = lognpdf(gc_grid, 0,green_dev_var);
green_EUR_dist = -100 .* gc_grid ./ theta_price(1) + green_lower_EUR;
green_wtp_1 = plot(green_EUR_dist,gc_dist);
title('Willingness to pay for green electricity')
xlabel('WTP (in EUR per month)');
ylabel('Density');
saveas(green_wtp_1,project_paths('OUT_FIGURES',['wtp_green_1_',file_suffix,'.pdf']));

X_EUR = 5;
coeff_X_EUR = - X_EUR .* theta_price(1) ./ 100;
share_green_larger_X = 1.0 - logncdf(coeff_X_EUR-theta_green(1), 0,theta_green(2));
share_green_positive = 1.0 - logncdf(0-theta_green(1), 0,theta_green(2));
fprintf('Share of consumers with positive willingness to pay for green electricity:\n %.4f\n', share_green_positive);
fprintf('Share of consumers with willingness to pay for green electricity of more than %d EUR:\n %.4f\n', X_EUR, share_green_larger_X);
   
if model_type==1 || model_type==2 || model_type==6 || model_type==7 || model_type==8
    fprintf('Implied magnitude of switching costs for mean income consumer: %4.0f EUR.\n',...
            sc_EUR(1))
elseif model_type==3 || model_type==4
    fprintf('Implied magnitude of switching costs for mean income young and senior consumer, respectively: %4.0f EUR and %4.0f EUR.\n',...
            sc_EUR(1), sc_EUR(2))
end

fprintf('Implied magnitude of incumbent preference for non-seniors (mean income consumer) : %4.2f EUR \n', ...
        inc_EUR(1))
fprintf('Implied magnitude of incumbent preference for seniors (mean income consumer) :%4.2f EUR \n', ...
        inc_EUR(2))    
fprintf('Implied magnitude of mean green preference for mean income consumer: %4.2f EUR \n', ...
        green_EUR(1))


% Compute implied PCW search costs in EUR.
% Keep in mind that we restrict PCW usage cost to be positive.
kappa_hat_young = exp([pcw_use_macro(:,2) pcw_use_macro(:,3) pcw_use_macro(:,5)] * theta_kappa(1:3)); 
if model_type==2 || model_type==4
    kappa_hat_senior = exp([pcw_use_macro(:,2) pcw_use_macro(:,3) pcw_use_macro(:,5)] * theta_kappa(1:3) + theta_kappa(4)); 
elseif model_type==1 || model_type==3 || model_type==5 || model_type==7 || model_type==8
    kappa_hat_senior = kappa_hat_young;
elseif model_type==6
    kappa_hat_young = zeros(T,1);
    kappa_hat_senior = kappa_hat_young;
end
kappa_hat_young_EUR = -100 * kappa_hat_young ./ theta_price(1);
kappa_hat_senior_EUR = -100 * kappa_hat_senior ./ theta_price(1);
% Compute average estimate of PCW search cost.
kappa_hat_avg_EUR = kappa_hat_young_EUR .* (1.0-0.189) +  kappa_hat_senior_EUR .* 0.189;
kappa_out_EUR = [kappa_hat_avg_EUR, kappa_hat_young_EUR, kappa_hat_senior_EUR];
% Print mean statistics on PCW usage costs.
fprintf('Mean (over time) PCW usage costs for average, young and senior consumer, respectively:\n%4.2f EUR\n%4.2f EUR\n%4.2f EUR\n',mean(kappa_out_EUR(:,1)), mean(kappa_out_EUR(:,2)),mean(kappa_out_EUR(:,3)));
% Write evolution of PCW entry cost over time to csv-file. 
dlmwrite(project_paths('OUT_ANALYSIS',['kappa_EUR_1_',file_suffix,'.csv']),kappa_out_EUR);

%% Extract some model statistics for specification comparisons.
% For first-stage GMM.
wtp_comp_label = {'Incumbent - non-senior', 'Incumbent - senior','Green electricity','PCW search costs - non-senior', 'PCW search cost - seniors', 'Switching cost - non-senior','Switching cost - senior'};
if model_type==4
    wtp_comp = [inc_EUR; green_EUR;mean(kappa_hat_young_EUR);mean(kappa_hat_senior_EUR);sc_EUR;];
else
    wtp_comp = [inc_EUR; green_EUR;mean(kappa_hat_young_EUR);mean(kappa_hat_senior_EUR);repmat(sc_EUR,2,1)];
end
wtp_comp_table = table(wtp_comp, 'RowNames',wtp_comp_label,'VariableNames', {'mean_wtp'});
writetable(wtp_comp_table, project_paths('OUT_ANALYSIS',['wtp_comp_1_',file_suffix,'.csv']),'writeRowNames',true);

%% Sanity check on predicted PCW usage and churn rates.
% Reformat observed churn rates.
churn_ts = reshape(churn_rates,J,T)';
churn_ts = churn_ts(:,1);
time_grid = linspace(1,T,T)';
% Export raw data for goodness of fit graph.
varNames_gf = {'month','pcwObs','pcwPred','churnObs','churnPred'};
gf_raw_data_1 = table(time_grid, pcw_use_macro(:,1), pred_PCWAgg1d_1, ...
    churn_ts, pred_churnAgg1d_1, 'VariableNames', varNames_gf);
writetable(gf_raw_data_1, project_paths('OUT_ANALYSIS',['gf_raw_data_1_',file_suffix,'.csv']));
% Profile out fixed effects from PCW usage.
PCW_pred_error = pred_PCWAgg1d_1 - pcw_use_macro(:,1);
% Define several dummy viarables.
yr_idx = [repmat(1,11,1); repmat(2,12,1); repmat(3,12,1); repmat(4,12,1); repmat(5,6,1)];
syr_idx = [repmat(1,5,1); kron(linspace(2,9,8)',ones(6,1))];
yr_dummies = dummyvar(yr_idx);
syr_dummies = dummyvar(syr_idx);
pcw_pred_err_reg = fitlm(yr_dummies(:,1),PCW_pred_error);
pcw_err_fitted = pcw_pred_err_reg.Fitted;
pred_PCWAgg1d_clean = pred_PCWAgg1d_1 - pcw_err_fitted; 
% Compute MA3 of both observed and predicted PCW usage.
MA_window = 3;
pred_pcw_MA = movmean(pred_PCWAgg1d_1,MA_window);
pred_churn_MA = movmean(pred_churnAgg1d_1,MA_window);
obs_pcw_MA = movmean(pcw_use_macro(:,1),MA_window);
obs_churn_MA = movmean(churn_ts,MA_window);
% Similarly for churn rates.
churn_pred_error = pred_churnAgg1d_1 - churn_ts;
% Define several dummy viarables.
churn_pred_err_reg = fitlm(yr_dummies(:,1),churn_pred_error);
churn_err_fitted = churn_pred_err_reg.Fitted;
pred_churnAgg1d_clean = pred_churnAgg1d_1 - churn_err_fitted; 
% Exploratory plot for PCW usage.
subplot(3,1,1)
plot(time_grid(1:end),pred_PCWAgg1d_1(1:end),time_grid(1:end),pcw_use_macro(1:end,1))
title('PCW usage observed vs predicted')
legend('Predicted','Observed')
subplot(3,1,2)
plot(time_grid(1:end),pred_PCWAgg1d_clean(1:end),time_grid(1:end),pcw_use_macro(1:end,1))
title('PCW usage observed vs predicted (year FE profiled out)')
legend('Predicted','Observed')
subplot(3,1,3)
plot(time_grid(1:end),pred_pcw_MA,time_grid(1:end),obs_pcw_MA);
legend('Predicted','Observed');
title('PCW usage observed vs predicted - MA(3)');
saveas(gcf,project_paths('OUT_FIGURES',['PCW_gf_1_',file_suffix,'.pdf']))

% Exploratory plot for churn rates.
subplot(3,1,1)
plot(time_grid(2:end),pred_churnAgg1d_1(1:end-1),time_grid(2:end),churn_ts(2:end,1))
title('Churn rates observed vs predicted')
legend('Predicted','Observed')
subplot(3,1,2)
plot(time_grid(2:end),pred_churnAgg1d_clean(1:end-1),time_grid(2:end),churn_ts(2:end,1))
title('Churn rates observed vs predicted (year FE profiled out)')
legend('Predicted','Observed ')
subplot(3,1,3)
plot(time_grid(2:end),pred_churn_MA(1:end-1),time_grid(2:end),obs_churn_MA(2:end));
legend('Predicted','Observed');
title('Churn rates observed vs predicted - MA(3)');
saveas(gcf,project_paths('OUT_FIGURES',['churn_gf_1_',file_suffix,'.pdf']))
% END OF EXPLORING GOF-STATISTICS FOR PCW USAGE AND CHURN RATES.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Various tables for illustrating model goodness of fit and plausibility.
% Create matrix of average valuations for different contracts and firms.
contract_id = linspace(1,J,J)';
firm_id = [1;1;2;2;3;4;4;5;5;6;7];
green_id = X_reg(1:J,3);
firm_values = firm_fe(1:J,:) * [theta_linear(1);theta_linear(3:end);0];
green_values =  green_id .* (theta_green(1) + green_dev_mean);
contract_values = firm_values + green_values;
contract_values_EUR = -100 .* contract_values ./ theta_price(1);
varNames_contract_values = {'firm','contract','green','coeff','eur'};
contract_values_out = table(firm_id, contract_id, green_id, contract_values, contract_values_EUR,'VariableNames', varNames_contract_values);
% Save average valuations of contracts after first-stage estimation.
writetable(contract_values_out, project_paths('OUT_ANALYSIS',['contract_values_1_',file_suffix,'.csv']),'writeRowNames',true);

% Create matrix of goodness of fit statistics for PCW usage and churn rates by year.
% Update only these variables for second stage.
gf_pcw_pred = pred_PCWAgg1d_1;
gf_churn_pred = pred_churnAgg1d_1;
gf_pcw_obs = pcw_use_macro(:,1);
gf_churn_obs = churn_ts;
% Matrix for churn goodness of fit.
%{
 Row 1: Observed mean PCW usage rate / churn rate
 Row 2: Predicted mean PCW usage rate / churn rate
 Row 3: Observed median PCW usage rate / churn rate
 Row 4: Predicted median PCW usage rate / churn rate
 Row 5: RMSE for PCW usage rate / churn rate
 Columns 1-5: for the five years in our sample
 Column 6: averaged over whole sample period
%}
gf_pcw_usage = zeros(5,6);
gf_churn_rates = zeros(5,6);
% Loop over years to fill all elements.
for yr=1:5
    yr_start = (yr-1)*12+1;
    yr_end = min(yr*12,T);
    % Fill elements for PCW usage statistics.
    gf_pcw_usage(1,yr) = mean(gf_pcw_obs(yr_start:yr_end));
    gf_pcw_usage(2,yr) = mean(gf_pcw_pred(yr_start:yr_end));
    gf_pcw_usage(3,yr) = median(gf_pcw_obs(yr_start:yr_end));
    gf_pcw_usage(4,yr) = median(gf_pcw_pred(yr_start:yr_end));
    % Compute RMSE for each year.
    gf_pcw_usage(5,yr) = sqrt(mean((gf_pcw_obs(yr_start:yr_end)-gf_pcw_pred(yr_start:yr_end)).^2));
    % Fill elements for churn rate statistics.
    gf_churn_rates(1,yr) = mean(gf_churn_obs(yr_start:yr_end));
    gf_churn_rates(2,yr) = mean(gf_churn_pred(yr_start:yr_end));
    gf_churn_rates(3,yr) = median(gf_churn_obs(yr_start:yr_end));
    gf_churn_rates(4,yr) = median(gf_churn_pred(yr_start:yr_end));
    % Compute RMSE for each year.
    gf_churn_rates(5,yr) = sqrt(mean((gf_churn_obs(yr_start:yr_end)-gf_churn_pred(yr_start:yr_end)).^2));
end
% Finally, write statistics for complete sample period.
gf_pcw_usage(1,end) = mean(gf_pcw_obs);
gf_pcw_usage(2,end) = mean(gf_pcw_pred);
gf_pcw_usage(3,end) = median(gf_pcw_obs);
gf_pcw_usage(4,end) = median(gf_pcw_pred);
gf_pcw_usage(5,end) = sqrt(mean((gf_pcw_obs-gf_pcw_pred).^2));
gf_churn_rates(1,end) = mean(gf_churn_obs);
gf_churn_rates(2,end) = mean(gf_churn_pred);
gf_churn_rates(3,end) = median(gf_churn_obs);
gf_churn_rates(4,end) = median(gf_churn_pred);
gf_churn_rates(5,end) = sqrt(mean((gf_churn_obs-gf_churn_pred).^2));
% Output matrices to csv-file for nice formatting in Python.
fprintf('Goodness of fit statistics by year for PCW usage:\n');
gf_pcw_usage
fprintf('Goodness of fit statistics by year for churn rates:\n');
gf_churn_rates
dlmwrite([output_path_log,['gf_pcw_1_',file_suffix,'.csv']],gf_pcw_usage,'delimiter',',');
dlmwrite([output_path_log,['gf_churn_1_',file_suffix,'.csv']],gf_churn_rates,'delimiter',',');

%% Call model file to compute CCP matrix and export to file.
run_estimation=0;
[~,~,~,~,~,~,~, ~,~, ~,~,~,~,~, pred_CCPiAvg_contract_1] = gmm_function(theta_opt);
run_estimation=1;
% Aggregate CCP output to firm level for each consumer.
pred_CCPiAvg_firm_1 = [ ...
    pred_CCPiAvg_contract_1(:,1,:) + pred_CCPiAvg_contract_1(:,2,:), ...
    pred_CCPiAvg_contract_1(:,3,:) + pred_CCPiAvg_contract_1(:,4,:), ...
    pred_CCPiAvg_contract_1(:,5,:) ...
    pred_CCPiAvg_contract_1(:,6,:) + pred_CCPiAvg_contract_1(:,7,:), ...
    pred_CCPiAvg_contract_1(:,8,:) + pred_CCPiAvg_contract_1(:,9,:), ...
    pred_CCPiAvg_contract_1(:,10,:) ...
    pred_CCPiAvg_contract_1(:,11,:)];
pred_CCPiAvg_firm_1([2,4,7,9],:,:) = [];

% Set list of simulated consumers to print for diagnostics.
consumer_list = [1;5;15;79;330];
fprintf('Example CCP matrices (averaged over time):\n');
pred_CCPiAvg_contract_1(:,:,consumer_list(1));
pred_CCPiAvg_contract_1(:,:,consumer_list(2));
pred_CCPiAvg_contract_1(:,:,consumer_list(3));
pred_CCPiAvg_contract_1(:,:,consumer_list(4));
pred_CCPiAvg_contract_1(:,:,consumer_list(5));
% Same on firm level.
pred_CCPiAvg_firm_1(:,:,consumer_list(1));
pred_CCPiAvg_firm_1(:,:,consumer_list(2));
pred_CCPiAvg_firm_1(:,:,consumer_list(3));
pred_CCPiAvg_firm_1(:,:,consumer_list(4));
pred_CCPiAvg_firm_1(:,:,consumer_list(5));
% Export CCPs for each consumer type averaged over time into
% three-dimensional arrays.
save([output_path_log,['ccp_contract_1_',file_suffix,'.mat']],'pred_CCPiAvg_contract_1');
save([output_path_log,['ccp_firm_1_',file_suffix,'.mat']],'pred_CCPiAvg_firm_1');

%% Compute standard errors and p-values.
% Calculate standard errors for first stage.
% Specify convergence criteria for contraction mapping routine.
% This can be commented out, to get SEs associated with a the convergence criteria that we use in the estimation. Since we only report SEs from the second stage below, this doesn't affect anything in the paper.
tol_delta_c = 1E-2;
max_iter_c = 90;
gmm_function = @(theta) model_myopic(theta, cdid, cdindex, market_shares, s_init_long, avg_churn_rates, mm_data_reshaped, pcw_use_micro, pcw_inst_micro, pcw_use_macro, adv_data, adv_data_mm, data_demo, X_reg, Z_matrix, draws_p_normalized, draws_eps,  preference_shocks,awareness_shocks,WMatrix, tol_delta_c,max_iter_c);

% Call function to compute standard errors.
if isequal(se_method,'simple')
    % Choose perturbation size that leads to well-behaved finite differences.
    stp_size = 1E-2;
    % Choose whether additive or multiplicative perturbation.
    ptb_mode = 'add';
    [varcov, SE] = calc_se(theta_opt, gmm_function, WMatrix, stp_size, ptb_mode); 
elseif isequal(se_method,'detail')
    % This part is not used in the final version.
    % Alternative to above simple computation of SE: compute SE for grid of
    % perturbations and multiplicative and additive perturbation types.
    % This grid may need adjustments.
    ptb_grid = [1e-1;1e-2;1e-3;1e-4;1e-5;1e-6];
    % Initialize containers for variance matrices.
    varcov_1_mult_all = zeros(length(theta_opt)+length(theta_linear),length(theta_opt)+length(theta_linear),length(ptb_grid));
    varcov_1_add_all = varcov_1_mult_all;    
    % Containers for comparing different SEs.
    std_errors = [];
    std_errors_mult =[];             
    std_errors_add  =[]; 
    % Loop over step sizes and perturbation types.
    for stp_index=1:length(ptb_grid)
        % Extract step size.
        stp_size = ptb_grid(stp_index);
        % Additive perturbation.
        [vcov_add, se_raw_add] = calc_se(theta_opt, gmm_function, WMatrix, stp_size, 'add');
        varcov_1_add_all(:,:,stp_index) = vcov_add;
        se_add = full(sqrt(diag(vcov_add)));
        std_errors_add = [std_errors_add;se_add'];
        % Multiplicative perturbation.
        [vcov_mult, se_raw_mult] = calc_se(theta_opt, gmm_function, WMatrix, stp_size, 'mult');
        varcov_1_mult_all(:,:,stp_index) = vcov_mult;
        se_mult = full(sqrt(diag(vcov_mult)));
        std_errors_mult = [std_errors_mult;se_mult'];
    end
    % Write "intermediate" element of SE sequence with additive perturbations into our standard container.
    vcov = varcov_1_mult_all(:,:,2);
    % Let's use a column vector here.
    SE = std_errors_mult(2,:)';
    std_errors = [std_errors;SE'];     
end
% Calculate t-statistics.
t_stats = [theta_linear; theta_opt] ./ SE;
% Calculate p-values.
p_values = 2.0 * (1.0 - normcdf(abs(t_stats)));    


%% For printing in results table: transform estimates for green coefficient into actual mean of distribution.   
% Compute mean of green coefficient distribution.
    green_mean_coeff = theta_green(1) + exp( theta_green(2).^2 ./2);
    % THIS DEPENDS ON MODEL SPECIFICATION.
    % Model 1: 2 - 12 (for green) & 1 - 13 (for incumbent-senior);
    % Model 4: 2 - 13 (for green) & 1 - 14 (for incumbent-senior);
    % Model 5: 2 - 12 (for green) & 1 - 13 (for incumbent-senior);
    % Model 6: 2 - 9 (for green) & 1 - 10 (for incumbent-senior);
    if model_type==1 || model_type==7 || model_type==8 % For model baseline.
        idx_green_lin = 2; % same across all mode specifications.
        idx_green_sd = 12; % varies depending on model specification.
        idx_inc_all = 1; % same across all mode specifications.
        idx_inc_senior = 13; % varies depending on model specification.
    elseif model_type==4
        idx_green_lin = 2; % same across all mode specifications.
        idx_green_sd = 12; % varies depending on model specification.
        idx_inc_all = 1; % same across all mode specifications.
        idx_inc_senior = 14; % varies depending on model specification.
    elseif model_type==5
        idx_green_lin = 2; % same across all mode specifications.
        idx_green_sd = 12; % varies depending on model specification.
        idx_inc_all = 1; % same across all mode specifications.
        idx_inc_senior = 13; % varies depending on model specification.
    elseif model_type==6
        idx_green_lin = 2; % same across all mode specifications.
        idx_green_sd = 12; % varies depending on model specification.
        idx_inc_all = 1; % same across all mode specifications.
        idx_inc_senior = 10; % varies depending on model specification.
    end
    % These computations do not depend on model specification.
    gprime_sigma = exp(theta_green(2)^2/2) * theta_green(2);
    delta_method_green = [0;1;zeros(idx_green_sd-3,1);gprime_sigma;zeros(length(SE)-idx_green_sd,1)]; 
    sd_green_mean =  sqrt(delta_method_green' * varcov * delta_method_green);
    gprime_inc_senior = 1;
    delta_method_inc_senior = [1;zeros(idx_inc_senior-2,1);1;zeros(length(SE)-idx_inc_senior,1)]; 
    sd_inc_senior = sqrt(delta_method_inc_senior' * varcov * delta_method_inc_senior);
    
    t_stat_green_mean = green_mean_coeff ./ sd_green_mean;
    p_value_green_mean = 2.0 * (1.0 - normcdf(abs(t_stat_green_mean)));   
    green_mean_out = [green_mean_coeff, sd_green_mean, p_value_green_mean,green_EUR];
    t_stat_inc_senior = theta_incumbent(2) ./ sd_inc_senior;
    p_value_inc_senior = 2.0 * (1.0 - normcdf(abs(t_stat_inc_senior)));   
    inc_senior_out = [theta_incumbent(2), sd_inc_senior, p_value_inc_senior,inc_EUR(2)];

%% Save estimation results to file.
theta_out = [theta_linear; theta_opt];
wtp_out = -100 .* theta_out ./ theta_price(1);

est_results = [theta_out, SE, p_values, wtp_out];

% This depends on model specification.
% Model 1: 2 (for green-sigma) - 13 (incumbent-seniors);
% Model 4: 2 (for green-sigma) - 14 (incumbent-seniors);
% Model 5: 2 (for green-sigma) - 13 (incumbent-seniors);
% Model 6: 2 (for green-sigma) - 10 (incumbent-seniors);
if model_type==1 || model_type==7 || model_type==8
    est_results(2,:) = green_mean_out;
    est_results(13,:) = inc_senior_out;
elseif model_type==4
    est_results(2,:) = green_mean_out;
    est_results(14,:) = inc_senior_out;
elseif model_type==5
    est_results(2,:) = green_mean_out;
    est_results(13,:) = inc_senior_out;
elseif model_type==6
    est_results(2,:) = green_mean_out;
    est_results(10,:) = inc_senior_out;
end

col_labels = {'Point Estimates', 'Standard Error','P-Values', 'Magnitude in EUR'};
if model_type==1 || model_type==7 || model_type==8
    row_labels = {...
                  'Incumbent (non-seniors)', ...
                  'Mean green coefficient', ...
                  'FE EDF', ...
                  'FE Eneco', ...
                  'FE ENINuon', ...
                  'FE Essent', ...
                  'FE Lampiris', ...
                  'PCW search', ...
                  'PCW search-Internet ', ...
                  'PCW search-Campaign', ...
                  'Mean price coefficient', ...
                  'Variance green coefficient', ...
                  'Incumbent (seniors)', ...
                  'Income-price interaction', ...
                  'Switching cost', ...
                  'Adv. constant', ...
                  'Adv. expenditure', ...
                  };
elseif model_type==2
    row_labels = {...
                  'Incumbent (non-seniors)', ...
                  'Mean green coefficient', ...
                  'FE EDF', ...
                  'FE Eneco', ...
                  'FE ENINuon', ...
                  'FE Essent', ...
                  'FE Lampiris', ...
                  'PCW search', ...
                  'PCW search-Internet ', ...
                  'PCW search-Campaign', ...
                  'PCW search-Senior', ...
                  'Mean price coefficient', ...
                  'Variance green coefficient', ...
                  'Incumbent (seniors)', ...
                  'Income-price interaction', ...
                  'Switching cost', ...
                  'Adv. constant', ...
                  'Adv. expenditure', ...
                  };
elseif model_type==3
    row_labels = {...
                  'Incumbent (non-seniors)', ...
                  'Mean green coefficient', ...
                  'FE EDF', ...
                  'FE Eneco', ...
                  'FE ENINuon', ...
                  'FE Essent', ...
                  'FE Lampiris', ...
                  'PCW search', ...
                  'PCW search-Internet ', ...
                  'PCW search-Campaign', ...
                  'Mean price coefficient', ...
                  'Variance green coefficient', ...
                  'Incumbent (seniors)', ...
                  'Income-price interaction', ...
                  'Switching cost', ...
                  'Switching cost (seniors)', ...
                  'Adv. constant', ...
                  'Adv. expenditure', ...
                  };
elseif model_type==4
    row_labels = {...
                  'Incumbent (non-seniors)', ...
                  'Mean green coefficient', ...
                  'FE EDF', ...
                  'FE Eneco', ...
                  'FE ENINuon', ...
                  'FE Essent', ...
                  'FE Lampiris', ...
                  'PCW search', ...
                  'PCW search-Internet ', ...
                  'PCW search-Campaign', ...
                  'PCW search-Senior', ...
                  'Mean price coefficient', ...
                  'Variance green coefficient', ...
                  'Incumbent (seniors)', ...
                  'Income-price interaction', ...
                  'Switching cost', ...
                  'Switching cost (seniors)', ...
                  'Adv. constant', ...
                  'Adv. expenditure', ...
                  };
elseif model_type==5
    row_labels = {...
                  'Incumbent (non-seniors)', ...
                  'Mean green coefficient', ...
                  'FE EDF', ...
                  'FE Eneco', ...
                  'FE ENINuon', ...
                  'FE Essent', ...
                  'FE Lampiris', ...
                  'PCW search', ...
                  'PCW search-Internet ', ...
                  'PCW search-Campaign', ...
                  'Mean price coefficient', ...
                  'Variance green coefficient', ...
                  'Incumbent (seniors)', ...
                  'Income-price interaction', ...
                  'Adv. constant', ...
                  'Adv. expenditure', ...
                  };
elseif model_type==6
    row_labels = {...
                  'Incumbent (non-seniors)', ...
                  'Mean green coefficient', ...
                  'FE EDF', ...
                  'FE Eneco', ...
                  'FE ENINuon', ...
                  'FE Essent', ...
                  'FE Lampiris', ...
                  'Mean price coefficient', ...
                  'Variance green coefficient', ...
                  'Incumbent (seniors)', ...
                  'Income-price interaction', ...
                  'Switching cost', ...
                  'Adv. constant', ...
                  'Adv. expenditure', ...
                  };
end             
% Write results to file.
% Save estimation results in raw format to be formatted in a nicer table using Python.
dlmwrite([output_path_log,['estresultsraw_1_',file_suffix,'.csv']],est_results);
fid = fopen([output_path_log,['estresultslegend_1_',file_suffix,'.csv']], 'w');
fprintf(fid, '%s,', row_labels{:});
fclose(fid);

%% Save full workspace as basis for counterfactuals.
save([output_path_log,['postestimation_workspace_1_',file_suffix,'.mat']]);
end
end

%% Run 2-step GMM estimation with efficient weighting matrix.
if efficient_GMM==1
    % Load relevant workspace file form first stage.
    load([output_path_log,['postestimation_workspace_1_',file_suffix,'.mat']]);
    % Set a few switches to make sure that relevant parts of model file code are run in the estimation.
    run_estimation = 1;
    SE_calc = 0;
    run_counterfactuals = 0;
    theta_opt_2_ars = theta_opt;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate efficient weighting matrix.
    WM_eff = calc_efficient_w(theta_opt,gmm_function);
    % Write relevant weighting matrix to global.
    WMatrix = WM_eff;
    WMEff_CN = cond(WMatrix);
    fprintf('Condition number of second-step weighting matrix is %e.\n',WMEff_CN);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Re-estimate with efficient weighting matrix.
    % Update gmm function with efficient weighting matrix.
    % Load convergence criteria for contraction mapping routine.
    load(project_paths('IN_MODEL_SPECS',['c_options_',file_suffix,'2']));
    
    gmm_function = @(theta) model_myopic(theta, cdid, cdindex, market_shares, s_init_long, avg_churn_rates, mm_data_reshaped, pcw_use_micro, pcw_inst_micro, pcw_use_macro, adv_data, adv_data_mm, data_demo, X_reg, Z_matrix, draws_p_normalized, draws_eps,  preference_shocks,awareness_shocks, WMatrix, tol_delta_c,max_iter_c);

    % (1) fminunc is gradient-based, therefore very efficient, but can also be
    % tricky if gradient does not behave well (likely to be the
    % case in our model).
    %[theta_opt, gmm_value] = fminsearch(gmm_function, theta_0,fminunc_options);
    % Run model file to update linear parameters.
    %gmm_function(theta_opt);
           
    % (2) ARS (accelerated random search): used in initial stages for finding rough list
    % of local minima under a broad search radius. Not used anymore.
    % Set contraction parameter.
    % Update search radius around estiamtes from first stage.
    % For debugging and exploration of parameter space: 
    LB = theta_opt - 0.5;
    UB = theta_opt + 0.5;
    c = 2;
    % Set convergence criterion.
    rho = 10^-6;
    % Here, the actual minimization routine is called.
    %theta_opt_2 = ars(gmm_function,LB,UB, c,rho,1000,6);
    % Transpose because ars.m outputs row vector.
    %theta_opt_2 = theta_opt_2';
    %theta_opt_2_ars = theta_opt_2;
        
    % (3) fminsearch: This is what we end up using.
    %% Load fminsearch options
    load(project_paths('IN_MODEL_SPECS',['fminsearch_options_',file_suffix,'2']));        
    format long
    % Run second stage estimation.
    [theta_opt_2, fval_opt_2, opt_flag_2, opt_out_2] = fminsearch(gmm_function, theta_opt_2_ars, fminsearch_options);
    format short
        
    % Save nonlinear parameters after second stage estimation.
    dlmwrite([output_path_log,['theta_opt_2_nl_',file_suffix,'.csv']],theta_opt_2,'delimiter',',');
    % Save linear parameters after second stage estimation.
    load(project_paths('OUT_ANALYSIS',['thetalinear_',file_suffix,'.mat']),'theta_linear');
    theta_lin_opt_2 = theta_linear;
    save(project_paths('OUT_ANALYSIS',['theta_lin_opt_2_',file_suffix]),'theta_lin_opt_2');
    % Save mean utilities after second stage estimation.
    load(project_paths('OUT_ANALYSIS',['expmvalold_',file_suffix,'.mat']),'expmvalold');
    expmval_opt_2 = expmvalold;
    save(project_paths('OUT_ANALYSIS',['expmvalold_opt_2_',file_suffix]),'expmval_opt_2');
    
    %% Some diagnostics in case of convergence issues.
    pshares_track = [market_shares];
    cSurplusPCW_track = [];
    cSurplusNoPCW_track = [];
    expmvalold_track = expmval_opt_2;
    theta_linear_track = theta_lin_opt_2;
    % Do another model evaluation at optimal parameter values to ensure convergence
    % and get predicted market shares, churn rates, and PCW usage.
    % This should not change anything.
    tol_delta_c = 1E-6;
    max_iter_c = 750;
        gmm_function = @(theta) model_myopic(theta, cdid, cdindex, market_shares, s_init_long, avg_churn_rates, mm_data_reshaped, pcw_use_micro, pcw_inst_micro, pcw_use_macro, adv_data, adv_data_mm, data_demo, X_reg, Z_matrix, draws_p_normalized, draws_eps,  preference_shocks,awareness_shocks,WMatrix, tol_delta_c,max_iter_c);
    % Do another model evaluation at optimal parameter values to ensure convergence
    % and get predicted market shares, churn rates, and PCW usage.
    % This should not change anything.
    [~,~,~,~,~,~,~, ...
        pred_shares1d_2, pred_churnAgg1d_2, ...
        ~,~, ...
        pred_PCWAgg1d_2,~,~, ...
        pred_CCPAvg_2] = gmm_function(theta_opt_2);
    % Check that linear parameters and mean utilities have not changed.
    load(project_paths('OUT_ANALYSIS',['thetalinear_',file_suffix,'.mat']),'theta_linear');
    load(project_paths('OUT_ANALYSIS',['expmvalold_',file_suffix,'.mat']),'expmvalold');
    conv_check_meanutil_2 = norm(expmval_opt_2-expmvalold);
    conv_check_thetalin_2 = norm(theta_lin_opt_2-theta_linear);
    % Track changes to various model statistics.
    pshares_track = [pshares_track, pred_shares1d_2];
    expmvalold_track = [expmvalold_track,expmvalold];
    theta_linear_track = [theta_linear_track, theta_linear];
    % Check changes across iterations.
    change_pshares = pshares_track(:,end) ./ pshares_track(:,end-1);
    change_expmval = expmvalold_track(:,end) ./ expmvalold_track(:,end-1);
    change_thetalinear = theta_linear_track(:,end) ./ theta_linear_track(:,end-1);

    % Do a final model evaluation at optimal parameter values to ensure convergence
    % and get predicted market shares, churn rates, and PCW usage.
    % This should not change anythign either.
    tol_delta_c = 1E-8;
    max_iter_c = 1000;
    gmm_function = @(theta) model_myopic(theta, cdid, cdindex, market_shares, s_init_long, avg_churn_rates, mm_data_reshaped, pcw_use_micro, pcw_inst_micro, pcw_use_macro, adv_data, adv_data_mm, data_demo, X_reg, Z_matrix, draws_p_normalized, draws_eps,  preference_shocks,awareness_shocks,WMatrix, tol_delta_c,max_iter_c);
    [~,~,~,~,~,~,~, ...
    pred_shares1d_2, pred_churnAgg1d_2, ...
    pred_cSurplusPCW_2,pred_cSurplusNoPCW_2, ...
    pred_PCWAgg1d_2,~,~, ...
    pred_CCPAvg_2] = gmm_function(theta_opt);
    
    % Check that linear parameters and mean utilities have not changed.
    load(project_paths('OUT_ANALYSIS',['thetalinear_',file_suffix,'.mat']),'theta_linear');
    load(project_paths('OUT_ANALYSIS',['expmvalold_',file_suffix,'.mat']),'expmvalold');
    conv_check_meanutil_2 = norm(expmval_opt_2-expmvalold);
    conv_check_thetalin_2 = norm(theta_lin_opt_2-theta_linear);

    % Track changes to various model statistics.
    pshares_track = [pshares_track, pred_shares1d_2];
    expmvalold_track = [expmvalold_track,expmvalold];
    theta_linear_track = [theta_linear_track, theta_linear];
    % Check changes across iterations.
    change_pshares = pshares_track(:,end) ./ pshares_track(:,end-1);
    change_expmval = expmvalold_track(:,end) ./ expmvalold_track(:,end-1);
    change_thetalinear = theta_linear_track(:,end) ./ theta_linear_track(:,end-1);

    %% Extract coefficients.
    % As for first-stage, this depends on the exact model specification that is run.
    if model_type==1 || model_type==7 || model_type==8
        theta_kappa = theta_opt_2(1:3);
        theta_price = [theta_opt_2(4); theta_opt_2(7)];
        theta_incumbent = [theta_linear(1); theta_linear(1) + theta_opt_2(6)];
        theta_green = [theta_linear(2); abs(theta_opt_2(5))];
        theta_sc = theta_opt_2(8);
        theta_adv = theta_opt_2(9:end);
    elseif model_type==2
        theta_kappa = theta_opt_2(1:4);
        theta_price = [theta_opt_2(5); theta_opt_2(8)];
        theta_incumbent = [theta_linear(1); theta_linear(1) + theta_opt_2(7)];
        theta_green = [theta_linear(2); abs(theta_opt_2(6))];
        theta_sc = theta_opt_2(9);
        theta_adv = theta_opt_2(10:end);
    elseif model_type==3
        theta_kappa = theta_opt_2(1:3);
        theta_price = [theta_opt_2(4); theta_opt_2(7)];
        theta_incumbent = [theta_linear(1); theta_linear(1) + theta_opt_2(6)];
        theta_green = [theta_linear(2); abs(theta_opt_2(5))];
        theta_sc = theta_opt_2(8:9);
        theta_adv = theta_opt_2(10:end);
    elseif model_type==4
        theta_kappa = theta_opt_2(1:4);
        theta_price = [theta_opt_2(5); theta_opt_2(8)];
        theta_incumbent = [theta_linear(1); theta_linear(1) + theta_opt_2(7)];
        theta_green = [theta_linear(2); abs(theta_opt_2(6))];
        theta_sc = theta_opt_2(9:10);
        theta_adv = theta_opt_2(11:end);
    elseif model_type==5
        theta_kappa = theta_opt_2(1:3);
        theta_price = [theta_opt_2(4); theta_opt_2(7)];
        theta_incumbent = [theta_linear(1); theta_linear(1) + theta_opt_2(6)];
        theta_green = [theta_linear(2); abs(theta_opt_2(5))];
        theta_sc = 0;
        theta_adv = theta_opt_2(8:end);
    elseif model_type==6
        theta_kappa = zeros(3,1);
        theta_price = [theta_opt_2(1); theta_opt_2(4)];
        theta_incumbent = [theta_linear(1); theta_linear(1) + theta_opt_2(3)];
        theta_green = [theta_linear(2); abs(theta_opt_2(2))];
        theta_sc = theta_opt_2(5);
        theta_adv = theta_opt_2(6:end);
    end
    %% Print results to screen.
        fprintf('Estimates for price coefficients for mean price and income-interaction coefficient: \n %6.4f \n %6.4f \n %6.4f \n', ...
                theta_price(1),theta_price(2))
        fprintf('Estimates for incumbent preference for young and senior: \n %6.4f \n %6.4f \n', ...
                theta_incumbent(1), theta_incumbent(2))
        fprintf('Estimates for green coefficients (mean and variance): \n ... !! CAREFUL WITH INTERPRETATION OF THIS IF PREFERENCE IS LOGNORMAL!! %6.4f \n %6.4f \n', ...
                theta_green(1), theta_green(2))
        if model_type==1 || model_type==2 || model_type==7 || model_type==8
            fprintf('Estimates for switching costs: \n %6.4f \n', ...
                theta_sc(1))
        elseif model_type==3 || model_type==4
            fprintf('Estimates for switching costs: \n %6.4f \n %6.4f \n', ...
                theta_sc(1), theta_sc(2))
        end
        fprintf('Estimates for advertising coefficients: \n %6.4f \n %6.4f \n\n ', ...
                theta_adv(1), theta_adv(2))

    % Calculate implied monetary magnitudes (for quantity model with log prices).
    % Willingness to pay for incumbent.
    inc_EUR = -100 .* theta_incumbent ./ repmat(theta_price(1),2,1);
    % Switching costs in EUR.
    sc_EUR = -100 * theta_sc ./ theta_price(1);

    % Willingness to pay for green electricity (mean).
    % Auxiliary computations for specification where deviation if modeled as
    % lognormal and total preference is given by "mean preference" (=lower
    % bound, theta_green(1)) + lognormal deviation with parameters 0 and theta_green(2).
    % Compute mean and variance of consumer-specific upward preference deviation.
    green_dev_mean = exp(0.5 * theta_green(2)^2);
    green_dev_var = exp( 2* theta_green(2)^2 - 1);
    % If green preference is specified as log-normal: Mean preference for green electricity.
    green_EUR = -100 ./ theta_price(1) .* (theta_green(1) + green_dev_mean);
    green_lower_EUR = -100 ./ theta_price(1) .* theta_green(1);
    % Only if green preference is specified as normal.
    % green_EUR = -100 .* theta_green(1) ./ theta_price(1);
    
    % Compute distribution of green electricity coefficient.
    % New version where green preference is modeled as log-normal distribution.
    gc_grid = linspace(0,15,100)';
    gc_dist = lognpdf(gc_grid, 0,green_dev_var);
    green_EUR_dist = -100 .* gc_grid ./ theta_price(1) + green_lower_EUR;
    subplot(1,1,1);
    green_wtp_2 = plot(green_EUR_dist,gc_dist);
    title('Willingness to pay for green electricity')
    xlabel('WTP (in EUR per month)');
    ylabel('Density');
    saveas(gcf,project_paths('OUT_FIGURES',['wtp_green_2_',file_suffix,'.pdf']));

    X_EUR = 5;
    coeff_X_EUR = - X_EUR .* theta_price(1) ./ 100;
    share_green_larger_X = 1.0 - logncdf(coeff_X_EUR-theta_green(1), 0,theta_green(2));
    share_green_positive = 1.0 - logncdf(0-theta_green(1), 0,theta_green(2));
    fprintf('Share of consumers with positive willingness to pay for green electricity:\n %.4f\n', share_green_positive);
    fprintf('Share of consumers with willingness to pay for green electricity of more than %d EUR:\n %.4f\n', X_EUR, share_green_larger_X);

    if model_type==1 || model_type==2 || model_type==7 || model_type==8
    fprintf('Implied magnitude of switching costs for mean income consumer: %4.0f EUR.\n',...
            sc_EUR(1))
    elseif model_type==3 || model_type==4
    fprintf('Implied magnitude of switching costs for mean income young and senior consumer, respectively: %4.0f EUR and %4.0f EUR.\n',...
            sc_EUR(1), sc_EUR(2))
    end
    fprintf('Implied magnitude of incumbent preference for non-seniors (mean income consumer) : %4.2f EUR \n', ...
            inc_EUR(1))
    fprintf('Implied magnitude of incumbent preference for seniors (mean income consumer) :%4.2f EUR \n', ...
            inc_EUR(2))    
    fprintf('Implied magnitude of mean green preference for mean income consumer: %4.2f EUR \n', ...
            green_EUR(1))

    % Compute implied PCW search costs in EUR.
    kappa_hat_young = exp([pcw_use_macro(:,2) pcw_use_macro(:,3) pcw_use_macro(:,5)] * theta_kappa(1:3)); 
    if model_type==1 || model_type==3 || model_type==5 || model_type==7 || model_type==8
        kappa_hat_senior = kappa_hat_young;
    elseif model_type==2 || model_type==4
        kappa_hat_senior = exp([pcw_use_macro(:,2) pcw_use_macro(:,3) pcw_use_macro(:,5)] * theta_kappa(1:3) + theta_kappa(4)); 
    end
    % Compute average statistics implied by PCW search cost estimates.
    kappa_hat_young_EUR = -100 * kappa_hat_young ./ theta_price(1);
    kappa_hat_senior_EUR = -100 * kappa_hat_senior ./ theta_price(1);
    kappa_hat_avg_EUR = kappa_hat_young_EUR .* (1.0-0.189) +  kappa_hat_senior_EUR .* 0.189;
    kappa_out_EUR = [kappa_hat_avg_EUR, kappa_hat_young_EUR, kappa_hat_senior_EUR];
    
    % Print mean statistics on PCW usage costs.
    fprintf('Mean (over time) PCW usage costs for average, young, senior, respectively:\n%6.4f EUR\n%6.4f EUR\n%6.4f EUR\n',mean(kappa_out_EUR(:,1)), mean(kappa_out_EUR(:,2)),mean(kappa_out_EUR(:,3)));
    % Sanity check on predicted PCW usage.
    time_grid = linspace(1,T,T)';
    plot(time_grid,pred_PCWAgg1d_2,time_grid,pcw_use_macro(:,1))
    legend('Predicted PCW usage','Observed PCW usage')
    % Write evolution of PCW entry cost over time to csv-file. 
    dlmwrite(project_paths('OUT_ANALYSIS',['kappa_EUR_2_',file_suffix,'.csv']),kappa_out_EUR);
    
    %% For efficient second-stage GMM.
    wtp_comp_label = {'Incumbent - non-senior', 'Incumbent - senior','Green electricity','PCW search costs - non-senior', 'PCW search cost - seniors', 'Switching cost - non-senior','Switching cost - senior'};
    if model_type==4
        wtp_comp = [inc_EUR; green_EUR;mean(kappa_hat_young_EUR);mean(kappa_hat_senior_EUR);sc_EUR];
    else
        wtp_comp = [inc_EUR; green_EUR;mean(kappa_hat_young_EUR);mean(kappa_hat_senior_EUR);repmat(sc_EUR,2,1);];
    end
    wtp_comp_table = table(wtp_comp, 'RowNames',wtp_comp_label,'VariableNames', {'mean_wtp'});
    writetable(wtp_comp_table, project_paths('OUT_ANALYSIS',['wtp_comp_2_',file_suffix,'.csv']),'writeRowNames',true);

    %% Sanity check on predicted PCW usage and churn rates after second stage GMM.
    % Reformat observed churn rates.
    churn_ts = reshape(churn_rates,J,T)';
    churn_ts = churn_ts(:,1);
    time_grid = linspace(1,T,T)';
    % Export raw data for goodness of fit graph.
    varNames_gf = {'month','pcwObs','pcwPred','churnObs','churnPred'};
    gf_raw_data_2 = table(time_grid, pcw_use_macro(:,1), pred_PCWAgg1d_2, ...
        churn_ts, pred_churnAgg1d_2, 'VariableNames', varNames_gf);
    writetable(gf_raw_data_2, project_paths('OUT_ANALYSIS',['gf_raw_data_2_',file_suffix,'.csv']));
    % Exploratory plot for goodness of fit.
    subplot(2,1,1)
    %plot(time_grid(2:end),pred_PCWAgg1d_1(2:end),time_grid(2:end),pcw_use_macro(1:end-1,1))
    plot(time_grid,pred_PCWAgg1d_2,time_grid,pcw_use_macro(:,1))
    legend('Predicted PCW usage','Observed PCW usage')
    subplot(2,1,2)
    %plot(time_grid(3:end),pred_churnAgg1d_1(2:end-1) ,time_grid(3:end),churn_ts(3:end))
    plot(time_grid(2:end),pred_churnAgg1d_2(1:end-1) ,time_grid(2:end),churn_ts(2:end))
    legend('Predicted churn rates','Observed churn rates')
    saveas(gcf,project_paths('OUT_FIGURES',['PCW_churn_gf_2_',file_suffix,'.pdf']))

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Various tables for illustrating model goodness of fit and plausibility.
    % Create matrix of average valuations for different contracts and firms.
    contract_id = linspace(1,J,J)';
    firm_id = [1;1;2;2;3;4;4;5;5;6;7];
    green_id = X_reg(1:J,3);
    firm_values = firm_fe(1:J,:) * [theta_linear(1);theta_linear(3:end);0];
    green_values =  green_id .* (theta_green(1) + green_dev_mean);
    contract_values = firm_values + green_values;
    contract_values_EUR = -100 .* contract_values ./ theta_price(1);
    varNames_contract_values = {'firm','contract','green','coeff','eur'};
    contract_values_out = table(firm_id, contract_id, green_id, contract_values, contract_values_EUR,'VariableNames', varNames_contract_values);
    % Save average valuations of contracts after first-stage estimation.
    writetable(contract_values_out, project_paths('OUT_ANALYSIS',['contract_values_2_',file_suffix,'.csv']),'writeRowNames',true);

    % Create matrix of goodness of fit statistics for PCW usage and churn rates by year.
    % Update only these variables for second stage.
    gf_pcw_pred = pred_PCWAgg1d_2;
    gf_churn_pred = pred_churnAgg1d_2;
    gf_pcw_obs = pcw_use_macro(:,1);
    gf_churn_obs = churn_ts;
    % Matrix for churn goodness of fit.
    %{
    Row 1: Observed mean PCW usage rate / churn rate
    Row 2: Predicted mean PCW usage rate / churn rate
    Row 3: Observed median PCW usage rate / churn rate
    Row 4: Predicted median PCW usage rate / churn rate
    Row 5: RMSE for PCW usage rate / churn rate
    Columns 1-5: for the five years in our sample
    Column 6: averaged over whole sample period
    %}
    gf_pcw_usage = zeros(5,6);
    gf_churn_rates = zeros(5,6);
    % Loop over years to fill all elements.
    for yr=1:5
        yr_start = (yr-1)*12+1;
        yr_end = min(yr*12,T);
        % Fill elements for PCW usage statistics.
        gf_pcw_usage(1,yr) = mean(gf_pcw_obs(yr_start:yr_end));
        gf_pcw_usage(2,yr) = mean(gf_pcw_pred(yr_start:yr_end));
        gf_pcw_usage(3,yr) = median(gf_pcw_obs(yr_start:yr_end));
        gf_pcw_usage(4,yr) = median(gf_pcw_pred(yr_start:yr_end));
        % Compute RMSE for each year.
        gf_pcw_usage(5,yr) = sqrt(mean((gf_pcw_obs(yr_start:yr_end)-gf_pcw_pred(yr_start:yr_end)).^2));
        % Fill elements for churn rate statistics.
        gf_churn_rates(1,yr) = mean(gf_churn_obs(yr_start:yr_end));
        gf_churn_rates(2,yr) = mean(gf_churn_pred(yr_start:yr_end));
        gf_churn_rates(3,yr) = median(gf_churn_obs(yr_start:yr_end));
        gf_churn_rates(4,yr) = median(gf_churn_pred(yr_start:yr_end));
        % Compute RMSE for each year.
        gf_churn_rates(5,yr) = sqrt(mean((gf_churn_obs(yr_start:yr_end)-gf_churn_pred(yr_start:yr_end)).^2));
    end
    % Finally, write statistics for complete sample period.
    gf_pcw_usage(1,end) = mean(gf_pcw_obs);
    gf_pcw_usage(2,end) = mean(gf_pcw_pred);
    gf_pcw_usage(3,end) = median(gf_pcw_obs);
    gf_pcw_usage(4,end) = median(gf_pcw_pred);
    gf_pcw_usage(5,end) = sqrt(mean((gf_pcw_obs-gf_pcw_pred).^2));
    gf_churn_rates(1,end) = mean(gf_churn_obs);
    gf_churn_rates(2,end) = mean(gf_churn_pred);
    gf_churn_rates(3,end) = median(gf_churn_obs);
    gf_churn_rates(4,end) = median(gf_churn_pred);
    gf_churn_rates(5,end) = sqrt(mean((gf_churn_obs-gf_churn_pred).^2));
    % Output matrices to csv-file for nice formatting in Python.
    fprintf('Goodness of fit statistics by year for PCW usage:\n');
    gf_pcw_usage
    fprintf('Goodness of fit statistics by year for churn rates:\n');
    gf_churn_rates
    dlmwrite([output_path_log,['gf_pcw_2_',file_suffix,'.csv']],gf_pcw_usage,'delimiter',',');
    dlmwrite([output_path_log,['gf_churn_2_',file_suffix,'.csv']],gf_churn_rates,'delimiter',',');

%% Call model file to compute CCP matrix and export to file (only for online appendix, and refere responses).
run_estimation=0;
[~,~,~,~,~,~,~, ~,~, ~,~,~,~,~, pred_CCPiAvg_contract_2] = gmm_function(theta_opt);
run_estimation=1;
% Aggregate CCP output to firm level for each consumer.
pred_CCPiAvg_firm_2 = [ ...
    pred_CCPiAvg_contract_2(:,1,:) + pred_CCPiAvg_contract_2(:,2,:), ...
    pred_CCPiAvg_contract_2(:,3,:) + pred_CCPiAvg_contract_2(:,4,:), ...
    pred_CCPiAvg_contract_2(:,5,:) ...
    pred_CCPiAvg_contract_2(:,6,:) + pred_CCPiAvg_contract_2(:,7,:), ...
    pred_CCPiAvg_contract_2(:,8,:) + pred_CCPiAvg_contract_2(:,9,:), ...
    pred_CCPiAvg_contract_2(:,10,:) ...
    pred_CCPiAvg_contract_2(:,11,:)];
pred_CCPiAvg_firm_2([2,4,7,9],:,:) = [];

% Set list of simulated consumers to print for diagnostics.
consumer_list = [1;5;15;79;330];
fprintf('Example CCP matrices (averaged over time):\n');
pred_CCPiAvg_contract_2(:,:,consumer_list(1));
pred_CCPiAvg_contract_2(:,:,consumer_list(2));
pred_CCPiAvg_contract_2(:,:,consumer_list(3));
pred_CCPiAvg_contract_2(:,:,consumer_list(4));
pred_CCPiAvg_contract_2(:,:,consumer_list(5));
% Same on firm level.
pred_CCPiAvg_firm_2(:,:,consumer_list(1));
pred_CCPiAvg_firm_2(:,:,consumer_list(2));
pred_CCPiAvg_firm_2(:,:,consumer_list(3));
pred_CCPiAvg_firm_2(:,:,consumer_list(4));
pred_CCPiAvg_firm_2(:,:,consumer_list(5));
% Export CCPs for each consumer type averaged over time into three-dimensional arrays.
save([output_path_log,['ccp_contract_2_',file_suffix,'.mat']],'pred_CCPiAvg_contract_2');
save([output_path_log,['ccp_firm_2_',file_suffix,'.mat']],'pred_CCPiAvg_firm_2');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate standard errors.
    % Specify convergence criteria for contraction mapping routine.
    load(project_paths('IN_MODEL_SPECS',['c_options_',file_suffix,'2']));
    gmm_function = @(theta) model_myopic(theta, cdid, cdindex, market_shares, s_init_long, avg_churn_rates, mm_data_reshaped, pcw_use_micro, pcw_inst_micro, pcw_use_macro, adv_data, adv_data_mm, data_demo, X_reg, Z_matrix, draws_p_normalized, draws_eps,  preference_shocks,awareness_shocks,WMatrix, tol_delta_c,max_iter_c);

    % Call function to compute standard errors.
    if isequal(se_method,'simple')
    % Set parameters that lead to wel-behaved finite differences.
    if model_type==6 % Since the model without PCW search costs leads to quite different parameter estimates, a smaller perturbation worked better here.
        stp_size = 1E-8;
    else
        stp_size = 1E-6;
    end
    % Choose whether to use additive or multiplicative perturbation.
    ptb_mode = 'add'
    [vcov_2, SE_2] = calc_se(theta_opt_2, gmm_function, WMatrix, stp_size, ptb_mode); 
    elseif isequal(se_method,'detail')
    %% This part is not used in the final version.
    % Alternative to above simple computation of SE: compute SE for gird of
    % perturbations and multiplicative and additive perturbation types.
    % This grid can be adjusted depending on specific data and model setup.
    ptb_grid = [1e-2;1e-4;1e-6];
    % Initialize containers for variance matrices.
    varcov_2_mult_all = zeros(length(theta_opt_2)+length(theta_linear),length(theta_opt_2)+length(theta_linear),length(ptb_grid));
    varcov_2_add_all = varcov_2_mult_all;    
    % Containers for comparing different SEs.
    std_errors_mult_2=[];             
    std_errors_add_2=[]; 
    % Loop over step sizes and perturbation types.
    for stp_index=1:length(ptb_grid)
        % Extract step size.
        stp_size = ptb_grid(stp_index);
        % Additive perturbation.
        [vcov_add, se_raw_add] = calc_se(theta_opt_2, gmm_function, WMatrix, stp_size, 'add');
        varcov_2_add_all(:,:,stp_index) = vcov_add;
        se_add = full(sqrt(diag(vcov_add)));
        std_errors_add_2 = [std_errors_add_2;se_add'];
        % Multiplicative perturbation.
        [vcov_mult, se_raw_mult] = calc_se(theta_opt_2, gmm_function, WMatrix, stp_size, 'mult');
        varcov_2_mult_all(:,:,stp_index) = vcov_mult;
        se_mult = full(sqrt(diag(vcov_mult)));
        std_errors_mult_2 = [std_errors_mult_2;se_mult'];
    end
    % Write "intermediate" element of SE sequence with additive perturbations into our standard container.
    vcov_2 = varcov_2_mult_all(:,:,2);
    % Let's use a column vector here.
    SE_2 = std_errors_mult_2(2,:)';
    % std_errors_2 = [std_errors_2;SE'];     
    end
    t_stats_2 = [theta_linear; theta_opt_2] ./ SE_2;
    % Calculate p-values.
    p_values_2 = 2.0 * (1.0 - normcdf(abs(t_stats_2)));    

    % Compute mean of green coefficient distribution.
    green_mean_coeff_2 = theta_green(1) + exp( theta_green(2).^2 ./2);
    % This depends on model specification.
    % Model 1: 2 - 12 (for green) & 1 - 13 (for incumbent-senior);
    % Model 4: 2 - 13 (for green) & 1 - 14 (for incumbent-senior);
    % Model 5: 2 - 12 (for green) & 1 - 13 (for incumbent-senior);
    % Model 6: 2 - 9 (for green) & 1 - 10 (for incumbent-senior);
    if model_type==1 || model_type==7 || model_type==8 % For model baseline.
        idx_green_lin = 2; % same across all mode specifications.
        idx_green_sd = 12; % varies depending on model specification.
        idx_inc_all = 1; % same across all mode specifications.
        idx_inc_senior = 13; % varies depending on model specification.
    elseif model_type==4
        idx_green_lin = 2; % same across all mode specifications.
        idx_green_sd = 12; % varies depending on model specification.
        idx_inc_all = 1; % same across all mode specifications.
        idx_inc_senior = 14; % varies depending on model specification.
    elseif model_type==5
        idx_green_lin = 2; % same across all mode specifications.
        idx_green_sd = 12; % varies depending on model specification.
        idx_inc_all = 1; % same across all mode specifications.
        idx_inc_senior = 13; % varies depending on model specification.
    elseif model_type==6
        idx_green_lin = 2; % same across all mode specifications.
        idx_green_sd = 12; % varies depending on model specification.
        idx_inc_all = 1; % same across all mode specifications.
        idx_inc_senior = 10; % varies depending on model specification.
    end
    % These computations do not depend on model specification.
    gprime_sigma = exp(theta_green(2)^2/2) * theta_green(2);
    delta_method_green = [0;1;zeros(idx_green_sd-3,1);gprime_sigma;zeros(length(SE)-idx_green_sd,1)]; 
    sd_green_mean =  sqrt(delta_method_green' * vcov_2 * delta_method_green);
    gprime_inc_senior = 1;
    delta_method_inc_senior = [1;zeros(idx_inc_senior-2,1);1;zeros(length(SE)-idx_inc_senior,1)]; 
    sd_inc_senior = sqrt(delta_method_inc_senior' * vcov_2 * delta_method_inc_senior);
    t_stat_green_mean = green_mean_coeff_2 ./ sd_green_mean;
    p_value_green_mean = 2.0 * (1.0 - normcdf(abs(t_stat_green_mean)));   
    green_mean_out = [green_mean_coeff_2, sd_green_mean, p_value_green_mean,green_EUR];

    t_stat_inc_senior = theta_incumbent(2) ./ sd_inc_senior;
    p_value_inc_senior = 2.0 * (1.0 - normcdf(abs(t_stat_inc_senior)));   
    inc_senior_out = [theta_incumbent(2), sd_inc_senior, p_value_inc_senior,inc_EUR(2)];

    
    
    %% Save estimation results to file.
    theta_out_2 = [theta_linear; theta_opt_2];
    wtp_out_2 = -100 .* theta_out_2 ./ theta_price(1);
    est_results_2 = [theta_out_2, SE_2, p_values_2, wtp_out_2];
    
    % This depends also on model specification.
    % Model 1: 2 (for green-sigma) - 13 (incumbent-seniors);
    % Model 4: 2 (for green-sigma) - 14 (incumbent-seniors);
    % Model 5: 2 (for green-sigma) - 13 (incumbent-seniors);
    % Model 6: 2 (for green-sigma) - 10 (incumbent-seniors);
    if model_type==1 || model_type==7 || model_type==8
        est_results_2(2,:) = green_mean_out;
        est_results_2(13,:) = inc_senior_out;
    elseif model_type==4
        est_results_2(2,:) = green_mean_out;
        est_results_2(14,:) = inc_senior_out;
    elseif model_type==5
        est_results_2(2,:) = green_mean_out;
        est_results_2(13,:) = inc_senior_out;
    elseif model_type==6
        est_results_2(2,:) = green_mean_out;
        est_results_2(10,:) = inc_senior_out;
    end
    col_labels = {'Point Estimates', 'Standard Error','P-Values', 'Magnitude in EUR'};
    % Define legends for parameter estimates.
    if model_type==1 || model_type==5 || model_type==7 || model_type==8
        row_labels = {...
                  'Incumbent (non-seniors)', ...
                  'Mean green coefficient', ...
                  'FE EDF', ...
                  'FE Eneco', ...
                  'FE ENINuon', ...
                  'FE Essent', ...
                  'FE Lampiris', ...
                  'PCW search', ...
                  'PCW search-Internet ', ...
                  'PCW search-Campaign', ...
                  'Mean price coefficient', ...
                  'Variance green coefficient', ...
                  'Incumbent (seniors)', ...
                  'Income-price interaction', ...
                  'Switching cost', ...
                  'Adv. constant', ...
                  'Adv. expenditure', ...
                  };
    elseif model_type==2
        row_labels = {...
                  'Incumbent (non-seniors)', ...
                  'Mean green coefficient', ...
                  'FE EDF', ...
                  'FE Eneco', ...
                  'FE ENINuon', ...
                  'FE Essent', ...
                  'FE Lampiris', ...
                  'PCW search', ...
                  'PCW search-Internet ', ...
                  'PCW search-Campaign', ...
                  'PCW search-Senior', ...
                  'Mean price coefficient', ...
                  'Variance green coefficient', ...
                  'Incumbent (seniors)', ...
                  'Income-price interaction', ...
                  'Switching cost', ...
                  'Adv. constant', ...
                  'Adv. expenditure', ...
                  };
    elseif model_type==3
        row_labels = {...
                  'Incumbent (non-seniors)', ...
                  'Mean green coefficient', ...
                  'FE EDF', ...
                  'FE Eneco', ...
                  'FE ENINuon', ...
                  'FE Essent', ...
                  'FE Lampiris', ...
                  'PCW search', ...
                  'PCW search-Internet ', ...
                  'PCW search-Campaign', ...
                  'Mean price coefficient', ...
                  'Variance green coefficient', ...
                  'Incumbent (seniors)', ...
                  'Income-price interaction', ...
                  'Switching cost', ...
                  'Switching cost (seniors)', ...
                  'Adv. constant', ...
                  'Adv. expenditure', ...
                  };
    elseif model_type==4
                row_labels = {...
                  'Incumbent (non-seniors)', ...
                  'Mean green coefficient', ...
                  'FE EDF', ...
                  'FE Eneco', ...
                  'FE ENINuon', ...
                  'FE Essent', ...
                  'FE Lampiris', ...
                  'PCW search', ...
                  'PCW search-Internet ', ...
                  'PCW search-Campaign', ...
                  'PCW search-Senior', ...
                  'Mean price coefficient', ...
                  'Variance green coefficient', ...
                  'Incumbent (seniors)', ...
                  'Income-price interaction', ...
                  'Switching cost', ...
                  'Switching cost (seniors)', ...
                  'Adv. constant', ...
                  'Adv. expenditure', ...
                  };
    elseif model_type==5
                row_labels = {...
                  'Incumbent (non-seniors)', ...
                  'Mean green coefficient', ...
                  'FE EDF', ...
                  'FE Eneco', ...
                  'FE ENINuon', ...
                  'FE Essent', ...
                  'FE Lampiris', ...
                  'PCW search', ...
                  'PCW search-Internet ', ...
                  'PCW search-Campaign', ...
                  'Mean price coefficient', ...
                  'Variance green coefficient', ...
                  'Incumbent (seniors)', ...
                  'Income-price interaction', ...
                  'Adv. constant', ...
                  'Adv. expenditure', ...
                  };
    elseif model_type==6
                row_labels = {...
                  'Incumbent (non-seniors)', ...
                  'Mean green coefficient', ...
                  'FE EDF', ...
                  'FE Eneco', ...
                  'FE ENINuon', ...
                  'FE Essent', ...
                  'FE Lampiris', ...
                  'Mean price coefficient', ...
                  'Variance green coefficient', ...
                  'Incumbent (seniors)', ...
                  'Income-price interaction', ...
                  'Switching cost', ...
                  'Adv. constant', ...
                  'Adv. expenditure', ...
                  };
    end
              
    % Save estimation results in raw format to be formatted in a nicer table using Python.
    dlmwrite([output_path_log,['estresultsraw_2_',file_suffix,'.csv']],est_results_2);
    fid = fopen([output_path_log,['estresultslegend_2_',file_suffix,'.csv']], 'w') ;
    fprintf(fid, '%s,', row_labels{:});
    fclose(fid) ;
    %% Save full workspace as basis for counterfactuals.
    save([output_path_log,['postestimation_workspace_2_',file_suffix,'.mat']]);
end
end % end loop over model types.
% Turn of log file.
diary OFF