"""
    Write several tables with model diagnostics information into tex-table.
    1. Average valuation of each contract
    2. Goodness-of-fit by Year for churn rates and PCW usage
    3. CCP matrices

    OUTPUT TABLES (one separate file for each combination of model_type and gmm_type, unless otherwise indicated):
    - gof_mm_pcw_FILESUFFIX.tex (we're probably not using those)
    - gof_mm_sw_ratio_FILESUFFIX.tex (we're probably not using those)
    
    - gof_churn_pcw_FILESUFFIX.tex

    - contract_val_FILESUFFIX.tex
    
    - ccp_contract_FILESUFFIX_CTYPE.tex (separate file for each simulated consumer types)
    - ccp_contract_old_sc_FILESUFFIX_CTYPE.tex (separate file for each simulated consumer types, only used for referee response)
    - ccp_firm_FILESUFFIX_CTYPE.tex (not used for paper, separate file for each simulated consumer types)
    - comp_contract_val_gmm_GMMTYPE (comparison across the four model specifications)
    - comp_est_stats_GMMTYPE (comparison across the four model specifications)

    OUTPUT FIGURES (one separate file for each combination of model_type and gmm_type):
    - cval_dist_FILESUFFIX.pdf


"""
import numpy as np
import pandas as pd
import scipy.io as sio
from bld.project_paths import project_paths_join
from matplotlib import pyplot as plt

# Select model types & GMM type to format.
model_list = [1, 4, 5, 6]
gmm_list = [
    [1, 2],
    [1, 2],
    [1],
    [1]]
model_print = r'Formatting estimation results for model {mt:d} and GMM stage {gmm:d}...\n'

# Initiate empty data frame for filling with contract valuations for comparison across model types.
model_names = ["Baseline", "Large", "No $\psi$", "No $\kappa$"]
# Data frame for comparing different contract valuations implied by different model specifications. Two data frames for two different GMM stages.
comp_contract_val_1 = pd.DataFrame(0, index=np.arange(11), columns=model_names)
comp_contract_val_2 = pd.DataFrame(0, index=np.arange(11), columns=model_names)
# For model comparison, we print 7 estimation statistics
comp_est_stats_1 = pd.DataFrame(0, index=np.arange(7), columns=model_names)
comp_est_stats_2 = pd.DataFrame(0, index=np.arange(7), columns=model_names)


for iter, model_type in enumerate(model_list):
    # print(iter)
    for gmm_type in gmm_list[iter]:
        print(model_print.format(mt=model_type, gmm=gmm_type))
        file_suffix = str(gmm_type) + '_mod' + str(model_type)
        file_suffix_modelonly = '1' + '_mod' + str(model_type)

        # Load observed and predicted micromoments for PCW usage and switching ratios as well as consumer-year weights.
        mm_weights = pd.read_csv(project_paths_join(
            'OUT_ANALYSIS', 'mm_data_weights_' + file_suffix_modelonly + '.csv'), header=None)
        pcw_active_obs = pd.read_csv(project_paths_join(
            'OUT_ANALYSIS', 'mm_data_active_' + file_suffix_modelonly + '.csv'), header=None)
        pcw_active_pred = pd.read_csv(project_paths_join(
            'OUT_ANALYSIS', 'mm_pred_active_' + file_suffix_modelonly + '.csv'), header=None)
        sw_ratio_obs = pd.read_csv(project_paths_join(
            'OUT_ANALYSIS', 'mm_data_swhet_' + file_suffix_modelonly + '.csv'), header=None)
        sw_ratio_pred = pd.read_csv(project_paths_join(
            'OUT_ANALYSIS', 'mm_pred_swhet_' + file_suffix_modelonly + '.csv'), header=None)
        # Compute prediction errors.
        pcw_active_err = pcw_active_obs-pcw_active_pred
        sw_ratio_err = sw_ratio_obs - sw_ratio_pred
        # Compute averages over year (unweighted).
        pcw_active_obs_avg_year = np.mean(pcw_active_obs, axis=1)
        pcw_active_pred_avg_year = np.mean(pcw_active_pred, axis=1)
        pcw_active_err_avg_year = np.mean(pcw_active_err, axis=1)
        sw_ratio_obs_avg_year = np.mean(sw_ratio_obs, axis=1)
        sw_ratio_pred_avg_year = np.mean(sw_ratio_pred, axis=1)
        sw_ratio_err_avg_year = np.mean(sw_ratio_err, axis=1)
        # Compute averages over consumer types (weighted).
        test = pcw_active_obs * mm_weights
        pcw_active_obs_avg_types = np.sum(pcw_active_obs * mm_weights, axis=0)
        pcw_active_pred_avg_types = np.sum(
            pcw_active_pred * mm_weights, axis=0)
        pcw_active_err_avg_types = np.sum(pcw_active_err * mm_weights, axis=0)
        sw_ratio_obs_avg_types = np.sum(sw_ratio_obs * mm_weights, axis=0)
        sw_ratio_pred_avg_types = np.sum(sw_ratio_pred * mm_weights, axis=0)
        sw_ratio_err_avg_types = np.sum(sw_ratio_err * mm_weights, axis=0)
        # Total average over both years and consumer types.
        pcw_active_obs_avg_total = np.mean(pcw_active_obs_avg_types, axis=0)
        pcw_active_pred_avg_total = np.mean(pcw_active_pred_avg_types, axis=0)
        pcw_active_err_avg_total = np.mean(pcw_active_err_avg_types, axis=0)
        sw_ratio_obs_avg_total = np.mean(sw_ratio_obs_avg_types, axis=0)
        sw_ratio_pred_avg_total = np.mean(sw_ratio_pred_avg_types, axis=0)
        sw_ratio_err_avg_total = np.mean(sw_ratio_err_avg_types, axis=0)

        # Define list of consumer types.
        ctypes = ['Inc. 1; non-senior',
                  'Inc. 2; non-senior',
                  'Inc. 3; non-senior',
                  'Inc. 4; non-senior',
                  'Inc. 1; senior',
                  'Inc. 2; senior',
                  'Inc. 3; senior',
                  'Inc. 4; senior'
                  ]

        # Goodness of fit by year: 2-panel table
        # Set different types of table rows with placeholders.
        table_row_header = '{var}  & & {val1} & {val2} & {val3} & {val4} & {val5} & {val6} \\tabularnewline \n'
        table_row = '{var} & {val1:.3f} & {val2:.3f} & {val3:.3f} & {val4:.3f} & {val5:.3f} & {val6:.3f} \\tabularnewline\n'
        table_row_type = """\multirow{{3}}{{*}}{{{var}}} & Obs. & {val1o:.4f} & {val2o:.4f} & {val3o:.4f} & {val4o:.4f} & {val5o:.4f} & {val6o:.4f} \\tabularnewline\n
            & Pred. & {val1p:.4f} & {val2p:.4f} & {val3p:.4f} & {val4p:.4f} & {val5p:.4f} & {val6p:.4f} \\tabularnewline\n
            & Err. & {val1e:.4f} & {val2e:.4f} & {val3e:.4f} & {val4e:.4f} & {val5e:.4f} & {val6e:.4f} \\tabularnewline\n """
        # First table: PCW usage (micro-level)
        # Second table: Switching ratios (micro-level)
        with open(project_paths_join('OUT_TABLES', 'gof_mm_pcw_' + file_suffix + '.tex'), 'w') as tex_file:
            # Top of table.
            tex_file.write('\\begin{tabular}{rr|ccccc|c}\n\\toprule\n')
            # Top panel for micro-level PCW usage rates.
            # Header row with names and dependent variables.
            tex_file.write(table_row_header.format(var='\\textbf{PCW usage}', val1='2012', val2='2013',
                                                   val3='2014', val4='2015',
                                                   val5='2016', val6='Avg (over years)'))
            tex_file.write('\\midrule \n')
            # Loop over consumer types (in row of array).
            # for index, row in gof_churn_yr.iterrows():
            for idx, row in pcw_active_obs.iterrows():
                # print(idx)
                # print(ctypes[idx])
                # print(row)
                tex_file.write(table_row_type.format(var=ctypes[idx],
                                                     val1o=pcw_active_obs[0][idx],
                                                     val2o=pcw_active_obs[1][idx],
                                                     val3o=pcw_active_obs[2][idx],
                                                     val4o=pcw_active_obs[3][idx],
                                                     val5o=pcw_active_obs[4][idx],
                                                     val6o=pcw_active_obs_avg_year[idx],
                                                     val1p=pcw_active_pred[0][idx],
                                                     val2p=pcw_active_pred[1][idx],
                                                     val3p=pcw_active_pred[2][idx],
                                                     val4p=pcw_active_pred[3][idx],
                                                     val5p=pcw_active_pred[4][idx],
                                                     val6p=pcw_active_pred_avg_year[idx],
                                                     val1e=pcw_active_err[0][idx],
                                                     val2e=pcw_active_err[1][idx],
                                                     val3e=pcw_active_err[2][idx],
                                                     val4e=pcw_active_err[3][idx],
                                                     val5e=pcw_active_err[4][idx],
                                                     val6e=pcw_active_err_avg_year[idx]
                                                     ))
                # Write separator line after each type.
                tex_file.write('\\midrule  \n')
            # Add row with average over types.
            tex_file.write(table_row_type.format(var='Avg. (over types)',
                                                 val1o=pcw_active_obs_avg_types[0],
                                                 val2o=pcw_active_obs_avg_types[1],
                                                 val3o=pcw_active_obs_avg_types[2],
                                                 val4o=pcw_active_obs_avg_types[3],
                                                 val5o=pcw_active_obs_avg_types[4],
                                                 val6o=pcw_active_obs_avg_total,
                                                 val1p=pcw_active_pred_avg_types[0],
                                                 val2p=pcw_active_pred_avg_types[1],
                                                 val3p=pcw_active_pred_avg_types[2],
                                                 val4p=pcw_active_pred_avg_types[3],
                                                 val5p=pcw_active_pred_avg_types[4],
                                                 val6p=pcw_active_pred_avg_total,
                                                 val1e=pcw_active_err_avg_types[0],
                                                 val2e=pcw_active_err_avg_types[1],
                                                 val3e=pcw_active_err_avg_types[2],
                                                 val4e=pcw_active_err_avg_types[3],
                                                 val5e=pcw_active_err_avg_types[4],
                                                 val6e=pcw_active_err_avg_total
                                                 ))
            # Write bottom of table and notes.
            tex_file.write('\\midrule \n')
            tex_file.write('\\midrule \n')
            tex_file.write(
                '\\multicolumn{8}{p{13.5cm}}{\\footnotesize{\\textit{Notes: The table displays the goodness of fit regarding PCW usage  on the micro (survey) level for each consumer type and year, as well as averaged in both dimensions.}}}\\tabularnewline\n')
            # Bottom of table.
            tex_file.write('\\bottomrule\n\\end{tabular}\n')
            # Second table  for switching ratio / heterogeneity.

        with open(project_paths_join('OUT_TABLES', 'gof_mm_sw_ratio_' + file_suffix + '.tex'), 'w') as tex_file:
            # Top of table.
            tex_file.write('\\begin{tabular}{rr|ccccc|c}\n\\toprule\n')
            # Top panel for micro-level PCW usage rates.
            # Header row with names and dependent variables.
            tex_file.write(table_row_header.format(var='\\textbf{Switching frequency}', val1='2012', val2='2013',
                                                   val3='2014', val4='2015',
                                                   val5='2016', val6='Avg (over years)'))
            tex_file.write('\\midrule \n')
            # Loop over consumer types (in row of array).
            # for index, row in gof_churn_yr.iterrows():
            for idx, row in sw_ratio_obs.iterrows():
                # Loop over consumer types (in row of array).
                tex_file.write(table_row_type.format(var=ctypes[idx],
                                                     val1o=sw_ratio_obs[0][idx],
                                                     val2o=sw_ratio_obs[1][idx],
                                                     val3o=sw_ratio_obs[2][idx],
                                                     val4o=sw_ratio_obs[3][idx],
                                                     val5o=sw_ratio_obs[4][idx],
                                                     val6o=sw_ratio_obs_avg_year[idx],
                                                     val1p=sw_ratio_pred[0][idx],
                                                     val2p=sw_ratio_pred[1][idx],
                                                     val3p=sw_ratio_pred[2][idx],
                                                     val4p=sw_ratio_pred[3][idx],
                                                     val5p=sw_ratio_pred[4][idx],
                                                     val6p=sw_ratio_pred_avg_year[idx],
                                                     val1e=sw_ratio_err[0][idx],
                                                     val2e=sw_ratio_err[1][idx],
                                                     val3e=sw_ratio_err[2][idx],
                                                     val4e=sw_ratio_err[3][idx],
                                                     val5e=sw_ratio_err[4][idx],
                                                     val6e=sw_ratio_err_avg_year[idx]
                                                     ))
                # Write separator line after each type.
                tex_file.write('\\midrule  \n')
            # Add row with average over types.
            tex_file.write(table_row_type.format(var='Avg. (over types)',
                                                 val1o=sw_ratio_obs_avg_types[0],
                                                 val2o=sw_ratio_obs_avg_types[1],
                                                 val3o=sw_ratio_obs_avg_types[2],
                                                 val4o=sw_ratio_obs_avg_types[3],
                                                 val5o=sw_ratio_obs_avg_types[4],
                                                 val6o=sw_ratio_obs_avg_total,
                                                 val1p=sw_ratio_pred_avg_types[0],
                                                 val2p=sw_ratio_pred_avg_types[1],
                                                 val3p=sw_ratio_pred_avg_types[2],
                                                 val4p=sw_ratio_pred_avg_types[3],
                                                 val5p=sw_ratio_pred_avg_types[4],
                                                 val6p=sw_ratio_pred_avg_total,
                                                 val1e=sw_ratio_err_avg_types[0],
                                                 val2e=sw_ratio_err_avg_types[1],
                                                 val3e=sw_ratio_err_avg_types[2],
                                                 val4e=sw_ratio_err_avg_types[3],
                                                 val5e=sw_ratio_err_avg_types[4],
                                                 val6e=sw_ratio_err_avg_total
                                                 ))
            # Write bottom of table and notes.
            tex_file.write('\\midrule  \n')
            tex_file.write('\\midrule  \n')
            tex_file.write(
                '\\multicolumn{8}{p{13.5cm}}{\\footnotesize{\\textit{Notes: The table displays the goodness of fit regarding relative switching frequency on a micro (survey) level for each consumer type and year, as well as averaged in both dimensions.}}}\\tabularnewline\n')
            # Bottom of table.
            tex_file.write('\\bottomrule\n\\end{tabular}\n')

            # Load estimation results for gross and net sample from csv-file.
            contract_values = pd.read_csv(project_paths_join(
                'OUT_ANALYSIS', 'contract_values_' + file_suffix + '.csv'))
            est_stats = pd.read_csv(project_paths_join(
                'OUT_ANALYSIS', 'wtp_comp_' + file_suffix + '.csv'))
            # Load row and column labels.
            contract_labels = ['ECS', 'ECS (g) ', 'EDF', 'EDF (g)', 'Essent (g)',
                               'ENI', 'ENI (g)', 'Eneco', 'Eneco (g)', 'Lampiris (g)', 'Other']
            firm_labels = ['ECS', 'EDF', 'Essent',
                           'ENI', 'Eneco', 'Lampiris', 'Other']

            # Compute valuations relative to incumbent conventional contract.
            coeff_ref = contract_values['coeff'][0]
            eur_ref = contract_values['eur'][0]
            contract_values['coeff'] = - coeff_ref + contract_values['coeff']
            contract_values['eur'] = - eur_ref + contract_values['eur']

            if gmm_type == 1:
                comp_contract_val_1.iloc[:, [iter]] = contract_values['eur']
                comp_est_stats_1.iloc[:, [iter]] = est_stats['mean_wtp']
            elif gmm_type == 2:
                comp_contract_val_2.iloc[:, [iter]] = contract_values['eur']
                comp_est_stats_2.iloc[:, [iter]] = est_stats['mean_wtp']

            # Table with average contract valuation.
            with open(project_paths_join('OUT_TABLES', 'contract_val_' + file_suffix + '.tex'), 'w') as tex_file:
                # Set different types of table rows with placeholders.
                table_row_header = '{val1} & {val2} & {val3} \\tabularnewline \n'
                table_row = '{name} & {val1:.2f} & {val2:.2f} \\tabularnewline\n'
                # Top of table.
                tex_file.write('\\begin{tabular}{r|cc}\n\\toprule\n')
                # Header row with names and dependent variables.
                tex_file.write(table_row_header.format(val1='Contract', val2='Valuation',
                                                       val3='WTP (in EUR/month)'))
                tex_file.write('\\midrule \n')
                # Write shares and HHI to table.
                for index, row in contract_values.iterrows():
                    tex_file.write(table_row.format(name=contract_labels[index],
                                                    val1=contract_values['coeff'][index],
                                                    val2=contract_values['eur'][index]))
                # Write bottom of table and notes.
                tex_file.write('\\midrule  \n')
                tex_file.write('\\multicolumn{3}{p{8cm}}{\\footnotesize{\\textit{Notes: The table displays the average valuation (across months and consumers) as measured by the sum of the relevant firm fixed effects and the mean preference for green electricity contracts.}}}\\tabularnewline\n')
                # Bottom of table.
                tex_file.write('\\bottomrule\n\\end{tabular}\n')

            # Goodness of fit variables.
            gof_vars = ['Mean (obs)', 'Mean (pred)', 'Median (obs)',
                        'Median (pred)', 'RMSE']
            gof_churn_yr = pd.read_csv(project_paths_join(
                'OUT_ANALYSIS', 'gf_churn_'+file_suffix + '.csv'), header=None)
            gof_pcw_yr = pd.read_csv(project_paths_join(
                'OUT_ANALYSIS', 'gf_pcw_'+file_suffix + '.csv'), header=None)

            # CCP matrices.
            ccp_contract = sio.loadmat(project_paths_join(
                'OUT_ANALYSIS', 'ccp_contract_' + file_suffix + '.mat'))
            ccp_firm = sio.loadmat(project_paths_join(
                'OUT_ANALYSIS', 'ccp_firm_' + file_suffix + '.mat'))
            ccp_contract_old_sc = sio.loadmat(project_paths_join(
                'IN_DATA', 'ccp_contract_old_sc_' + file_suffix + '.mat'))

            # Goodness of fit by year: 2-panel table
            # Top panel: Churn rates
            # Bottom panel: PCW usage
            # Table with average contract valuation.
            with open(project_paths_join('OUT_TABLES', 'gof_churn_pcw_' + file_suffix + '.tex'), 'w') as tex_file:
                # Set different types of table rows with placeholders.
                table_row_header = '{var} & {val1} & {val2} & {val3} & {val4} & {val5} & {val6} \\tabularnewline \n'
                table_row = '{var} & {val1:.3f} & {val2:.3f} & {val3:.3f} & {val4:.3f} & {val5:.3f} & {val6:.3f} \\tabularnewline\n'

                # Top of table.
                tex_file.write('\\begin{tabular}{r|ccccc|c}\n\\toprule\n')
                # Top panel for chunr rates.
                # Header row with names and dependent variables.
                tex_file.write(table_row_header.format(var='\\textbf{Churn rates}', val1='2012', val2='2013',
                                                       val3='2014', val4='2015',
                                                       val5='2016', val6='Average'))
                tex_file.write('\\midrule \n')
                # Write shares and HHI to table.
                for index, row in gof_churn_yr.iterrows():
                    # print(gof_vars[index])
                    if gof_vars[index] != 'RMSE':
                        tex_file.write(table_row.format(var=gof_vars[index],
                                                        val1=gof_churn_yr[0][index],
                                                        val2=gof_churn_yr[1][index],
                                                        val3=gof_churn_yr[2][index],
                                                        val4=gof_churn_yr[3][index],
                                                        val5=gof_churn_yr[4][index],
                                                        val6=gof_churn_yr[5][index]))
                # Write bottom of table and notes.
                tex_file.write('\\midrule  \n')
                tex_file.write('\\midrule  \n')
                # Bottom panel for PCW usage.
                # Header row with names and dependent variables.
                tex_file.write(table_row_header.format(var='\\textbf{PCW usage}', val1='2012', val2='2013',
                                                       val3='2014', val4='2015',
                                                       val5='2016', val6='Average'))
                tex_file.write('\\midrule \n')
                for index, row in gof_pcw_yr.iterrows():
                    if gof_vars[index] != 'RMSE':
                        tex_file.write(table_row.format(var=gof_vars[index],
                                                        val1=gof_pcw_yr[0][index],
                                                        val2=gof_pcw_yr[1][index],
                                                        val3=gof_pcw_yr[2][index],
                                                        val4=gof_pcw_yr[3][index],
                                                        val5=gof_pcw_yr[4][index],
                                                        val6=gof_pcw_yr[5][index]))
                # Write bottom of table and notes.
                tex_file.write('\\midrule  \n')
                tex_file.write('\\midrule  \n')
                tex_file.write(
                    '\\multicolumn{7}{p{12.5cm}}{\\footnotesize{\\textit{Notes: The table displays the goodness of fit of our baseline model regarding industry-level churn rates and PCW usage behavior.}}}\\tabularnewline\n')
                # Bottom of table.
                tex_file.write('\\bottomrule\n\\end{tabular}\n')

            # CCP matrices.
            ccp_contract_dict = sio.loadmat(project_paths_join(
                'OUT_ANALYSIS', 'ccp_contract_' + file_suffix + '.mat'))
            ccp_contract_old_sc_dict = sio.loadmat(project_paths_join(
                'IN_DATA', 'ccp_contract_old_sc_' + file_suffix + '.mat'))
            ccp_firm_dict = sio.loadmat(project_paths_join(
                'OUT_ANALYSIS', 'ccp_firm_' + file_suffix + '.mat'))

            # Relevant CCP matrices for representative consumer types.
            if gmm_type == 1:
                ccp_contract = ccp_contract_dict['pred_CCPiAvg_contract_1']
                ccp_contract_old_sc = ccp_contract_old_sc_dict['pred_CCPiAvg_contract_1']
                ccp_firm = ccp_firm_dict['pred_CCPiAvg_firm_1']
            elif gmm_type == 2:
                ccp_contract = ccp_contract_dict['pred_CCPiAvg_contract_2']
                ccp_contract_old_sc = ccp_contract_old_sc_dict['pred_CCPiAvg_contract_2']
                ccp_firm = ccp_firm_dict['pred_CCPiAvg_firm_2']

            # Indexes for four representative consumer types (based on MATLAB index, recall Python is 0-based); always include consumer 0 to have some dependency constant in waf-script.
            # Young consumer with low income: 0
            # Young consumer with mean income: 128
            # Young consumer with high income: 374
            # Senior consumer with low income: 200
            # Senior consumer with mean income: 278
            # Senior consumer with high income: 99
            consumer_idx = [0, 128, 374, 200, 278, 99]

            # Set different types of table rows with placeholders.
            table_row_header_contract = ' & {val1} & {val2} & {val3} & {val4} & {val5} & {val6} & {val7} & {val8} & {val9} & {val10} & {val11} \\tabularnewline \n'
            table_row_contract = '{var} & {val1:.4f} & {val2:.4f} & {val3:.4f} & {val4:.4f} & {val5:.4f} & {val6:.4f} & {val7:.4f} & {val8:.4f} & {val9:.4f} & {val10:.4f} & {val11:.4f} \\tabularnewline\n'
            table_row_header_firm = ' & {val1} & {val2} & {val3} & {val4} & {val5} & {val6} & {val7} \\tabularnewline \n'
            table_row_firm = '{var} & {val1:.4f} & {val2:.4f} & {val3:.4f} & {val4:.4f} & {val5:.4f} & {val6:.4f} & {val7:.4f} \\tabularnewline\n'

            # CCP matrices for different consumer types
            # Loop over consumer types.
            for c_type in consumer_idx:

                # Double-check this indeed results in two-dimensional array
                ccp_contract_i = ccp_contract[:, :, c_type]
                ccp_contract_i_old_sc = ccp_contract_old_sc[:, :, c_type]
                ccp_firm_i = ccp_firm[:, :, c_type]
                # Table with contract-level CCP matrix.
                with open(project_paths_join('OUT_TABLES', 'ccp_contract_' + str(c_type) + '_' + file_suffix + '.tex'), 'w') as tex_file:
                    # Top of table.
                    tex_file.write(
                        '\\begin{tabular}{r|ccccccccccc}\n\\toprule\n')
                    # Top panel for chunr rates.
                    # Header row with names and dependent variables.
                    tex_file.write(table_row_header_contract.format(val1=contract_labels[0], val2=contract_labels[1], val3=contract_labels[2],
                                                                    val4=contract_labels[3], val5=contract_labels[
                                                                        4], val6=contract_labels[5],
                                                                    val7=contract_labels[6], val8=contract_labels[
                                                                        7], val9=contract_labels[8],
                                                                    val10=contract_labels[9], val11=contract_labels[10]))
                    tex_file.write('\\midrule \n')
                    # Write shares and HHI to table.
                    # for index, row in ccp_contract_i.iterrows():
                    # for count, value in enumerate(values)
                    for idx, row in enumerate(ccp_contract_i):
                        # print(idx)
                        # print(row)
                        tex_file.write(table_row_contract.format(var=contract_labels[idx],
                                                                 val1=ccp_contract_i[idx][0],
                                                                 val2=ccp_contract_i[idx][1],
                                                                 val3=ccp_contract_i[idx][2],
                                                                 val4=ccp_contract_i[idx][3],
                                                                 val5=ccp_contract_i[idx][4],
                                                                 val6=ccp_contract_i[idx][5],
                                                                 val7=ccp_contract_i[idx][6],
                                                                 val8=ccp_contract_i[idx][7],
                                                                 val9=ccp_contract_i[idx][8],
                                                                 val10=ccp_contract_i[idx][9],
                                                                 val11=ccp_contract_i[idx][10]))
                    # Write bottom of table and notes.
                    tex_file.write('\\midrule  \n')
                    tex_file.write('\\midrule  \n')
                    tex_file.write(
                        '\\multicolumn{12}{p{19cm}}{\\footnotesize{\\textit{Notes: The table displays the predicted conditional choice probabilities (on the contract-level) averaged over time for one representative consumer type.}}}\\tabularnewline\n')
                    # Bottom of table.
                    tex_file.write('\\bottomrule\n\\end{tabular}\n')

                # Only for referee replies: Table with contract-level CCP matrix when switching costs are firm-specific.
                with open(project_paths_join('OUT_TABLES', 'ccp_contract_old_sc_' + str(c_type) + '_' + file_suffix + '.tex'), 'w') as tex_file:
                    # Top of table.
                    tex_file.write(
                        '\\begin{tabular}{r|ccccccccccc}\n\\toprule\n')
                    # Top panel for chunr rates.
                    # Header row with names and dependent variables.
                    tex_file.write(table_row_header_contract.format(val1=contract_labels[0], val2=contract_labels[1], val3=contract_labels[2],
                                                                    val4=contract_labels[3], val5=contract_labels[
                                                                        4], val6=contract_labels[5],
                                                                    val7=contract_labels[6], val8=contract_labels[
                                                                        7], val9=contract_labels[8],
                                                                    val10=contract_labels[9], val11=contract_labels[10]))
                    tex_file.write('\\midrule \n')
                    # Write shares and HHI to table.
                    # for index, row in ccp_contract_i.iterrows():
                    # for count, value in enumerate(values)
                    for idx, row in enumerate(ccp_contract_i_old_sc):
                        # print(idx)
                        # print(row)
                        tex_file.write(table_row_contract.format(var=contract_labels[idx],
                                                                 val1=ccp_contract_i_old_sc[idx][0],
                                                                 val2=ccp_contract_i_old_sc[idx][1],
                                                                 val3=ccp_contract_i_old_sc[idx][2],
                                                                 val4=ccp_contract_i_old_sc[idx][3],
                                                                 val5=ccp_contract_i_old_sc[idx][4],
                                                                 val6=ccp_contract_i_old_sc[idx][5],
                                                                 val7=ccp_contract_i_old_sc[idx][6],
                                                                 val8=ccp_contract_i_old_sc[idx][7],
                                                                 val9=ccp_contract_i_old_sc[idx][8],
                                                                 val10=ccp_contract_i_old_sc[idx][9],
                                                                 val11=ccp_contract_i_old_sc[idx][10]))
                    # Write bottom of table and notes.
                    tex_file.write('\\midrule  \n')
                    tex_file.write('\\midrule  \n')
                    tex_file.write(
                        '\\multicolumn{12}{p{19cm}}{\\footnotesize{\\textit{Notes: The table displays the predicted conditional choice probabilities (on the contract-level) averaged over time for one representative consumer type when switching costs are firm-specific (as in our initial version) and not contract-specific (as in the current version).}}}\\tabularnewline\n')
                    # Bottom of table.
                    tex_file.write('\\bottomrule\n\\end{tabular}\n')

                with open(project_paths_join('OUT_TABLES', 'ccp_firm_' + str(c_type) + '_' + file_suffix + '.tex'), 'w') as tex_file:

                    # Top of table.
                    tex_file.write('\\begin{tabular}{r|ccccccc}\n\\toprule\n')
                    # Top panel for chunr rates.
                    # Header row with names and dependent variables.
                    tex_file.write(table_row_header_firm.format(val1=firm_labels[0], val2=firm_labels[1], val3=firm_labels[2], val4=firm_labels[3], val5=firm_labels[4], val6=firm_labels[5],
                                                                val7=firm_labels[6]))
                    tex_file.write('\\midrule \n')

                    # Loop over contracts/firms.
                    for idx, row in enumerate(ccp_firm_i):
                        # print(idx)
                        # print(row)
                        tex_file.write(table_row_firm.format(var=firm_labels[idx],
                                                             val1=ccp_firm_i[idx][0],
                                                             val2=ccp_firm_i[idx][1],
                                                             val3=ccp_firm_i[idx][2],
                                                             val4=ccp_firm_i[idx][3],
                                                             val5=ccp_firm_i[idx][4],
                                                             val6=ccp_firm_i[idx][5],
                                                             val7=ccp_firm_i[idx][6]))
                    # Write bottom of table and notes.
                    tex_file.write('\\midrule  \n')
                    tex_file.write('\\midrule  \n')
                    tex_file.write(
                        '\\multicolumn{8}{p{11cm}}{\\footnotesize{\\textit{Notes: The table displays the predicted conditional choice probabilities (on the firm-level) averaged over time for one representative consumer type.}}}\\tabularnewline\n')
                    # Bottom of table.
                    tex_file.write('\\bottomrule\n\\end{tabular}\n')

        # contract_labels = ['ECS', 'ECS (g) ', 'EDF', 'EDF (g)', 'Essent (g)',
#                    'ENI', 'ENI (g)', 'Eneco', 'Eneco (g)', 'Lampiris (g)', 'Other']
# firm_labels = ['ECS', 'EDF', 'Essent', 'ENI', 'Eneco', 'Lampiris', 'Other']

        # Load distribution of contract valuations.
        cval_dist_dict = sio.loadmat(project_paths_join(
            'OUT_ANALYSIS', 'contract_val_dist_' + file_suffix + '.mat'))
        cval_dist = cval_dist_dict['v_i_all']
        # Extract dimensional of data.
        ns_cons = cval_dist.shape[2]
        T = cval_dist.shape[0]
        J = cval_dist.shape[1]
        # Move loop dimension to first to allow zip to work smoothly.
        cval_dist = np.moveaxis(cval_dist, 1, 0)

        # Create histograms.
        # Subplot 2 x 5 for distribution of each carrier.
        plt.clf()
        fig, axarr = plt.subplots(nrows=2, ncols=5, figsize=(11.5, 8))
        fig.subplots_adjust(hspace=0.5)
        fig.suptitle('Estimated Distribution of contract valuation for model ' +
                     str(model_type) + ' and GMM stage ' + str(gmm_type), fontsize=14)
        fig.subplots_adjust(hspace=0.5)

        # Extract mimum and maximum valuations for scaling x axis.
        cval_max = np.amax(cval_dist)
        cval_min = np.amin(cval_dist)

        for ax, data, name in zip(axarr.flatten(), cval_dist, contract_labels):
            # print(ax)
            # print(data.flatten().shape)
            # print(name)
            data_reshaped = data.flatten()
            ax.hist(data_reshaped, color='green', align='mid',
                    bins=25,  range=[cval_min, cval_max], lw=0)
            ax.grid(False)
            ax.axes.get_xaxis().set_visible(False)
            ax.axes.get_yaxis().set_visible(False)
            ax.set(title=name)
            ax.axvline(data_reshaped.mean(), color='k',
                       linestyle='dashed', linewidth=1)
            min_ylim, max_ylim = ax.get_ylim()

            ax.text(data_reshaped.mean()+1.0, max_ylim*0.85,
                    'Mean: {:.2f}'.format(data_reshaped.mean()))
        # Save the figure.
        plt.savefig(project_paths_join(
            'OUT_FIGURES', 'cval_dist_'+file_suffix + '.pdf'), bbox_inches='tight')


# Format table that compares contract valuations implied by different model specifications.
# Set different types of table rows with placeholders.
table_row_header_comp = ' {varname}  & {val1} & {val2} & {val3} & {val4} \\tabularnewline \n'
table_row_comp = '{name} & {val1: .2f} & {val2: .2f} & {val3: .2f} & {val4: .2f} \\tabularnewline\n'

for gmm_type in gmm_list[0]:
    # Extract relevant data.
    if gmm_type == 1:
        comp_data = comp_contract_val_1
    elif gmm_type == 2:
        comp_data = comp_contract_val_2

    # Table with average contract valuation.
    with open(project_paths_join('OUT_TABLES', 'comp_contract_val_gmm_' + str(gmm_type) + '.tex'), 'w') as tex_file:
        # Top of table.
        tex_file.write('\\begin{tabular}{r|cccc}\n\\toprule\n')
        # Header row with names and dependent variables.
        tex_file.write(table_row_header_comp.format(
            varname='Contract', val1=model_names[0], val2=model_names[1], val3=model_names[2], val4=model_names[3]))
        tex_file.write('\\midrule \n')
        # Write shares and HHI to table.
        for index, row in comp_data.iterrows():
            # print(index)
            if index < J-1:  # do not print outside good valuation.
                tex_file.write(table_row_comp.format(name=contract_labels[index],
                                                     val1=comp_data.iloc[index, 0],
                                                     val2=comp_data.iloc[index, 1],
                                                     val3=comp_data.iloc[index, 2],
                                                     val4=comp_data.iloc[index, 3]))
        # Write bottom of table and notes.
        tex_file.write('\\midrule  \n')
        tex_file.write('\\multicolumn{5}{p{8cm}}{\\footnotesize{\\textit{Notes: The table compares the average contract valuation (across months and consumers) as measured by the sum of the relevant firm fixed effects and the mean preference for green electricity contracts obtained from different model specifications. Valuations are measured in EUR per month relative to the conventional contract of the incumbent firm.}}}\\tabularnewline\n')
        # Bottom of table.
        tex_file.write('\\bottomrule\n\\end{tabular}\n')

    # Format comparison of important estimation statistics/WTPs from the four different model.
    # Extract relevant data.
    if gmm_type == 1:
        comp_est_stats = comp_est_stats_1
    elif gmm_type == 2:
        comp_est_stats = comp_est_stats_2

    comp_est_stats['Row'] = est_stats['Row']
    comp_est_stats = comp_est_stats.set_index('Row')
    # Table with average contract valuation.
    with open(project_paths_join('OUT_TABLES', 'comp_est_stats_' + str(gmm_type) + '.tex'), 'w') as tex_file:
        # Top of table.
        tex_file.write('\\begin{tabular}{r|cccc}\n\\toprule\n')
        # Header row with names and dependent variables.
        tex_file.write(table_row_header_comp.format(
            varname='Average WTP', val1=model_names[0], val2=model_names[1], val3=model_names[2], val4=model_names[3]))
        tex_file.write('\\midrule \n')
        # Write shares and HHI to table.
        for index, row in comp_est_stats.iterrows():
            #print(index)
            tex_file.write(table_row_comp.format(name=index,
                                                 val1=comp_est_stats['Baseline'][index],
                                                 val2=comp_est_stats['Large'][index],
                                                 val3=comp_est_stats['No $\psi$'][index],
                                                 val4=comp_est_stats['No $\kappa$'][index]))
        # Write bottom of table and notes.
        tex_file.write('\\midrule \n')
        tex_file.write('\\multicolumn{5}{p{12cm}}{\\footnotesize{\\textit{Notes: The table compares the average WTP (across months and consumers) for various product characteristics and consumer types obtained from different model specifications. Valuations are measured in EUR per month relative to the outside good.}}}\\tabularnewline\n')
        # Bottom of table.
        tex_file.write('\\bottomrule\n\\end{tabular}\n')


# Create final versions of results for paper.
# Set different types of table rows with placeholders.
table_row_header_comp = ' {varname}  & {val1} & {val2} & {val3} \\tabularnewline \n'
table_row_comp = '{name} & {val1: .2f} & {val2: .2f} & {val3: .2f} \\tabularnewline\n'
# Initiate empty data frame for filling with contract valuations for comparison across model types.
model_names = ["Full model", "No switching costs", "No PCW search costs"]

comp_data['Baseline'] = comp_contract_val_2['Baseline']
comp_data['No $\psi$'] = comp_contract_val_1['No $\psi$']
comp_data['No $\kappa$'] = comp_contract_val_1['No $\kappa$']

# Final table with average contract valuation.
with open(project_paths_join('OUT_TABLES', 'comp_contract_val_final' + '.tex'), 'w') as tex_file:
    # Top of table.
    tex_file.write('\\begin{tabular}{r|ccc}\n\\toprule\n')
    # Header row with names and dependent variables.
    tex_file.write(table_row_header_comp.format(
        varname='Contract', val1=model_names[0], val2=model_names[1], val3=model_names[2]))
    tex_file.write('\\midrule \n')
    # Write shares and HHI to table.
    for index, row in comp_data.iterrows():
        # print(index)
        if index < J-1:  # do not print outside good valuation.
            tex_file.write(table_row_comp.format(name=contract_labels[index],
                                                 val1=comp_data.iloc[index, 0],
                                                 val2=comp_data.iloc[index, 2],
                                                 val3=comp_data.iloc[index, 3]))
    # Write bottom of table and notes.
    tex_file.write('\\midrule  \n')
    tex_file.write('\\multicolumn{4}{p{12cm}}{\\footnotesize{\\textit{Notes: The table compares the average contract valuation (across months and consumers) as measured by the sum of the relevant firm fixed effects and the mean preference for green electricity contracts obtained from different model specifications. Valuations are measured in EUR per month relative to the conventional contract of the incumbent firm.}}}\\tabularnewline\n')
    # Bottom of table.
    tex_file.write('\\bottomrule\n\\end{tabular}\n')


# Format comparison of important estimation statistics/WTPs from the four different model.
# Extract relevant data.
comp_est_stats_1 = comp_est_stats_1.set_index('Row')
comp_est_stats_2 = comp_est_stats_2.set_index('Row')
comp_est_stats = comp_est_stats_2
#comp_est_stats['Row'] = est_stats['Row']
#comp_est_stats = comp_est_stats.set_index('Row')
comp_est_stats['Baseline'] = comp_est_stats_2['Baseline']
comp_est_stats['No $\psi$'] = comp_est_stats_1['No $\psi$']
comp_est_stats['No $\kappa$'] = comp_est_stats_1['No $\kappa$']

# Table with average contract valuation.
with open(project_paths_join('OUT_TABLES', 'comp_est_stats_final' + '.tex'), 'w') as tex_file:
    # Top of table.
    tex_file.write('\\begin{tabular}{r|ccc}\n\\toprule\n')
    # Header row with names and dependent variables.
    tex_file.write(table_row_header_comp.format(
        varname='Average WTP', val1=model_names[0], val2=model_names[1], val3=model_names[2]))
    tex_file.write('\\midrule \n')
    # Write shares and HHI to table.
    for index, row in comp_est_stats.iterrows():
        # print(index)
        tex_file.write(table_row_comp.format(name=index,
                                             val1=comp_est_stats['Baseline'][index],
                                             val2=comp_est_stats['No $\psi$'][index],
                                             val3=comp_est_stats['No $\kappa$'][index]))
    # Write bottom of table and notes.
    tex_file.write('\\midrule \n')
    tex_file.write('\\multicolumn{4}{p{16cm}}{\\footnotesize{\\textit{Notes: The table compares the average WTP (across months and consumers) for various product characteristics and consumer types obtained from different model specifications.}}}\\tabularnewline\n')
    # Bottom of table.
    tex_file.write('\\bottomrule\n\\end{tabular}\n')
