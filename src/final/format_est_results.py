"""
    Write the results of the estimation as table code into  a .tex-file
    for tables with complete estimation results for Appendix for all model types and GMM stages.

    Output tables (separate file for each combination of model_list and gmm_list):
    - out/tables/estresults_FILESUFFIX.tex

    Most importantly, this file formats Table 2 in the main text, as well as Figure 1 in the main text.


"""
import csv
import numpy as np
import pandas as pd
from bld.project_paths import project_paths_join
pd.options.mode.chained_assignment = None  # default='warn'


# Select model types & GMM type to format.
model_list = [1, 4, 5, 6, 7, 8]
gmm_list = [
    [1, 2],
    [1, 2],
    [1],
    [1],
    [1, 2],
    [1, 2]]
model_print = r'Formatting estimation results for model {mt:d} and GMM stage {gmm:d}...\n'

for iter, model_type in enumerate(model_list):
    for gmm_type in gmm_list[iter]:
        print(model_print.format(mt=model_type, gmm=gmm_type))
        # Define file suffix for current model specification.
        file_suffix = str(gmm_type) + '_mod' + str(model_type)

        # Load estimation results from csv-file.
        legend = pd.read_csv(project_paths_join(
            'OUT_ANALYSIS', 'estresultslegend_' + file_suffix + '.csv'), header=None, index_col=False)
        legend.drop(legend.columns[[-1]], axis=1, inplace=True)
        # Only for debugging, remove once output from estimation file is corrected.
        #legend.drop(legend.columns[[0]], axis=1, inplace=True)
        legend = legend.transpose().reset_index()

        eres = pd.read_csv(project_paths_join(
            'OUT_ANALYSIS', 'estresultsraw_' + file_suffix + '.csv'), header=None, index_col=False)
        eres.columns = ['beta', 'sigma', 'pval', 'WTP']
        # Drop linear firm fixed effects parameters.
        # eres.drop([2,3,4,5,6],inplace=True)
        # eres.reset_index()
        eres['var_name'] = legend[0]

        # Add significance stars.
        eres['sig_stars'] = ''
        eres.loc[eres.pval > 0.1, 'sig_stars'] = ''
        eres.loc[(eres.pval <= 0.1) & (eres.pval > 0.05), 'sig_stars'] = '*'
        eres.loc[(eres.pval <= 0.05) & (eres.pval > 0.01), 'sig_stars'] = '**'
        eres.loc[eres.pval <= 0.01, 'sig_stars'] = '***'


        # Copy objects for main text table.
        if model_type==1 and gmm_type==2:
            legend_1 = legend.copy(deep=True)
            eres_1 = eres.copy(deep=True)
        elif model_type==4 and gmm_type==2:
            legend_4 = legend.copy(deep=True)
            eres_4 = eres.copy(deep=True)
        elif model_type==7 and gmm_type==2:
            legend_7 = legend.copy(deep=True)
            eres_7 = eres.copy(deep=True)
        # Add significance stars.
        eres['sig_stars'] = ''
        eres.loc[eres.pval > 0.1, 'sig_stars'] = ''
        eres.loc[(eres.pval <= 0.1) & (eres.pval > 0.05), 'sig_stars'] = '*'
        eres.loc[(eres.pval <= 0.05) & (eres.pval > 0.01), 'sig_stars'] = '**'
        eres.loc[eres.pval <= 0.01, 'sig_stars'] = '***'


        # Set a different types of table row with placeholders.
        table_row_var = """\multirow{{2}}{{*}}{{{var}}} & {val1:.4f}{str1} & {val2: .2f} \\tabularnewline\n  & ({se1:.4f}) & \\tabularnewline\n"""
        #table_row_var_2 = """\multirow{{2}}{{*}}{{{var}}} &  & {val2:.4f}{str2} & {val3:.4f}{str3}  \\tabularnewline\n  &  & ({se2:.4f}) & ({se3:.4f}) \\tabularnewline\n"""
        #table_row_var_3 = """\multirow{{2}}{{*}}{{{var}}} &  &  & {val3:.4f}{str3}  \\tabularnewline\n  &  &  & ({se3:.4f}) \\tabularnewline\n"""
        #table_row_int = '{stat} & {val1:.4g} & {val2:.4g} & {val3:.4g}  \\tabularnewline\n'
        #table_row_float = '{stat} & {val1:.3f} & {val2:.3f} & {val3:.3f}  \\tabularnewline\n'
        table_row_strings = '{stat} & {val1} & {val2} \\tabularnewline \n'
        #table_row_test_1 = '{test} & {val1:.3f} & {val2:.3f} & {val3:.3f}  \\tabularnewline\n'
        #table_row_test_2 = '{test} &  & {val2:.3f} & {val3:.3f}  \\tabularnewline\n'
        #table_row_test_3 = '{test} &  & & {val3:.3f} \\tabularnewline\n'
        table_row_price = """\multirow{{2}}{{*}}{{{var}}} & {val1:.4f}{str1} & - \\tabularnewline\n  & ({se1:.4f}) & \\tabularnewline\n"""


        # Switch order of coefficients and drop meaningless WTPs from column.
        if model_type == 1 or model_type == 7 or model_type == 8:
            row_order = [10, 13, 0, 12, 1, 11, 14, 7, 8, 9, 15, 16]
            info_order = [7, 8, 9, 15, 16]
        elif model_type == 4:
            row_order = [11, 14, 0, 13, 1, 12, 15, 16, 7, 8, 9, 10, 17, 18]
            info_order = [7, 8, 9, 10, 17, 18]
        elif model_type == 5:
            row_order = [10, 13, 0, 12, 1, 11, 7, 8, 9, 14, 15]
            info_order = [7, 8, 9, 14, 15]
        elif model_type == 6:
            row_order = [7, 10, 0, 9, 1, 8, 11, 12, 13]
            info_order = [12, 13]

        # print('Order of original coefficients in table:')
        # print(row_order)
        # Write the results to a LaTeX table.
        with open(project_paths_join('OUT_TABLES', 'estresults_' + file_suffix + '.tex'), 'w') as tex_file:

            # Top of table.
            tex_file.write('\\begin{tabular}{lcc}\n\\toprule\n')
            # Header row with names and dependent variables.
        #     tex_file.write(table_row_strings.format(stat='', val1='Main Model',
        #                                            val2=''))
            tex_file.write(table_row_strings.format(stat='', val1='Coefficients',
                                                    val2='WTP in EUR'))
        #     tex_file.write(table_row_strings.format(stat='', val1=ols_results['reg_names'][0],
        #                                            val2=ols_results['reg_names'][1],
        #                                            val3=ols_results['reg_names'][2]))

            tex_file.write('\\midrule')
        # legend.at[0, row_idx]
            # Write coefficients to table.
            for i, var in enumerate(row_order):
                # row_idx = row_order[i]
                row_idx = var
               # print(row_idx)
                # if row_idx in(3,4,5,7,9,11,12):
                if model_type == 1 or model_type == 7 or model_type == 8:
                    if row_idx in (10, 13, 11, 7, 8, 9, 15, 16):
                        tex_file.write(table_row_price.format(var=eres['var_name'][row_idx], val1=eres['beta'][row_idx],
                                                              str1=eres['sig_stars'][row_idx],
                                                              se1=eres['sigma'][row_idx]))
                    else:
                        tex_file.write(table_row_var.format(var=eres['var_name'][row_idx], val1=eres['beta'][row_idx],
                                                            val2=eres['WTP'][row_idx],
                                                            str1=eres['sig_stars'][row_idx],
                                                            se1=eres['sigma'][row_idx]))
                elif model_type == 4:
                    if row_idx in (11, 14, 12, 7, 8, 9, 10, 17, 18):
                        tex_file.write(table_row_price.format(var=eres['var_name'][row_idx], val1=eres['beta'][row_idx],
                                                              str1=eres['sig_stars'][row_idx],
                                                              se1=eres['sigma'][row_idx]))
                    else:
                        tex_file.write(table_row_var.format(var=eres['var_name'][row_idx], val1=eres['beta'][row_idx],
                                                            val2=eres['WTP'][row_idx],
                                                            str1=eres['sig_stars'][row_idx],
                                                            se1=eres['sigma'][row_idx]))
                elif model_type == 5:
                    if row_idx in (10, 13, 11, 7, 8, 9, 14, 15):
                        tex_file.write(table_row_price.format(var=eres['var_name'][row_idx], val1=eres['beta'][row_idx],
                                                              str1=eres['sig_stars'][row_idx],
                                                              se1=eres['sigma'][row_idx]))
                    else:
                        tex_file.write(table_row_var.format(var=eres['var_name'][row_idx], val1=eres['beta'][row_idx],
                                                            val2=eres['WTP'][row_idx],
                                                            str1=eres['sig_stars'][row_idx],
                                                            se1=eres['sigma'][row_idx]))
                elif model_type == 6:
                    if row_idx in (7, 10, 8, 12, 13):
                        tex_file.write(table_row_price.format(var=eres['var_name'][row_idx], val1=eres['beta'][row_idx],
                                                              str1=eres['sig_stars'][row_idx],
                                                              se1=eres['sigma'][row_idx]))
                    else:
                        tex_file.write(table_row_var.format(var=eres['var_name'][row_idx], val1=eres['beta'][row_idx],
                                                            val2=eres['WTP'][row_idx],
                                                            str1=eres['sig_stars'][row_idx],
                                                            se1=eres['sigma'][row_idx]))

            tex_file.write('\\midrule\n')
            if gmm_type == 1:
                tex_file.write('\\multicolumn{3}{p{11cm}}{\\footnotesize{\\textit{Notes: Results from estimating the demand model using GMM with block-diagonal 2SLS weighting matrix. Standard errors in parentheses. *,**,*** denote significance at the 10\%, 5\% and 1\%-level respectively. - denotes non-interpretable willingness-to-pay.}}}\\tabularnewline\n')
            elif gmm_type == 2:
                tex_file.write('\\multicolumn{3}{p{11cm}}{\\footnotesize{\\textit{Notes: Results from estimating the demand model using efficient 2-step GMM  weighting matrix. Standard errors in parentheses. *,**,*** denote significance at the 10\%, 5\% and 1\%-level respectively. - denotes non-interpretable willingness-to-pay.}}}\\tabularnewline\n')
            # Bottom of table.
            tex_file.write('\\bottomrule\n\\end{tabular}\n')

        # Write separate table for PCW search costs and advertising coefficients only.
        if model_type != 6:
            # Update variable legend.
            if model_type == 1 or model_type == 7 or model_type == 8:
                eres['var_name'][7] = 'PCW search cost - constant'
                eres['var_name'][8] = 'PCW search cost - internet penetration rate'
                eres['var_name'][9] = 'PCW search cost - campaign'
                # eres['var_name'][10] = 'PCW search cost - senior'
                eres['var_name'][15] = 'Awareness process - constant'
                eres['var_name'][16] = 'Awareness process - adv. expenditure'
            elif model_type == 4:
                eres['var_name'][7] = 'PCW search cost - constant'
                eres['var_name'][8] = 'PCW search cost - internet penetration rate'
                eres['var_name'][9] = 'PCW search cost - campaign'
                eres['var_name'][10] = 'PCW search cost - senior'
                eres['var_name'][17] = 'Awareness process - constant'
                eres['var_name'][18] = 'Awareness process - adv. expenditure'
            elif model_type == 5:
                eres['var_name'][7] = 'PCW search cost - constant'
                eres['var_name'][8] = 'PCW search cost - internet penetration rate'
                eres['var_name'][9] = 'PCW search cost - campaign'
                eres['var_name'][14] = 'Awareness process - constant'
                eres['var_name'][15] = 'Awareness process - adv. expenditure'

            # Define new table row types.
            table_row_var_info = """\multirow{{2}}{{*}}{{{var}}} & {val1:.4f}{str1} \\tabularnewline\n  & ({se1:.4f}) \\tabularnewline\n"""
            # table_row_header_info = """\multirow{{2}}{{*}}{{{var}}} & {str1} \\tabularnewline\n """
            table_row_header_info = '{stat} & {val1}  \\tabularnewline \n'

            # # Write the results to a LaTeX table.
            # with open(project_paths_join('OUT_TABLES', 'estresultsinfoparameters_' + file_suffix + '.tex'), 'w') as tex_file:
            #     # Top of table.
            #     tex_file.write('\\begin{tabular}{lcc}\n\\toprule\n')
            #     # Header row with names and dependent variables.
            #     tex_file.write(table_row_header_info.format(
            #         stat='', val1='Coefficients'))
            #     tex_file.write('\\midrule')

            #     # Write coefficients to table.
            #     for i, row_idx in enumerate(info_order):
            #         tex_file.write(table_row_var_info.format(
            #             var=eres['var_name'][row_idx], val1=eres['beta'][row_idx], str1=eres['sig_stars'][row_idx], se1=eres['sigma'][row_idx]))
            #         # Midrule between two groups of parameters.
            #         if row_idx == 9 & model_type == 1:
            #             tex_file.write('\\midrule\n')
            #         elif row_idx == 9 & model_type == 5:
            #             tex_file.write('\\midrule\n')
            #         elif row_idx == 10 & model_type == 4:
            #             tex_file.write('\\midrule\n')

            #     # write bottom of table.
            #     tex_file.write('\\midrule\n')
            #     if gmm_type == 1:
            #         tex_file.write('\\multicolumn{2}{p{9.5cm}}{\\footnotesize{\\textit{Notes: Results for parameter for advertising-awareness process and PCW search costs from estimating the demand model using GMM with block-diagonal 2SLS weighting matrix. Standard errors in parentheses. *,**,*** denote significance at the 10\%, 5\% and 1\%-level respectively.}}}\\tabularnewline\n')
            #     elif gmm_type == 2:
            #         tex_file.write('\\multicolumn{2}{p{9.5cm}}{\\footnotesize{\\textit{Notes: Results for parameter for advertising-awareness process and PCW search costs from estimating the demand model using efficient two-step GMM. Standard errors in parentheses. *,**,*** denote significance at the 10\%, 5\% and 1\%-level respectively.}}}\\tabularnewline\n')
            #     # Bottom of table.
            #     tex_file.write('\\bottomrule\n\\end{tabular}\n')


# For referee responses, add comparison table for different logit smoother values.
model_print_comp = r'Formatting model comparison of model {mt1:d} (baseline) and model {mt2:d} (larger logit smoother)...\n'

print(model_print_comp.format(mt1=1, mt2=8))
# Define file suffix for current model specification.
# Baseline model.
file_suffix_base = str(1) + '_mod' + str(1)
# Baseline modle with larger logit smoother value.
file_suffix_lsv = str(1) + '_mod' + str(8)

# Load estimation results from csv-file.
legend = pd.read_csv(project_paths_join(
    'OUT_ANALYSIS', 'estresultslegend_' + file_suffix + '.csv'), header=None, index_col=False)
legend.drop(legend.columns[[-1]], axis=1, inplace=True)
# Only for debugging, remove once output from estimation file is corrected.
#legend.drop(legend.columns[[0]], axis=1, inplace=True)
legend = legend.transpose().reset_index()

eres_base = pd.read_csv(project_paths_join(
    'OUT_ANALYSIS', 'estresultsraw_' + file_suffix_base + '.csv'), header=None, index_col=False)
eres_base.columns = ['beta', 'sigma', 'pval', 'WTP']
eres_lsv = pd.read_csv(project_paths_join(
    'OUT_ANALYSIS', 'estresultsraw_' + file_suffix_lsv + '.csv'), header=None, index_col=False)
eres_lsv.columns = ['beta', 'sigma', 'pval', 'WTP']
eres_base['var_name'] = legend[0]
eres_lsv['var_name'] = legend[0]

# Set a different types of table row with placeholders.
table_row_var_comp = """\multirow{{2}}{{*}}{{{var}}} & {val1:.4f}{str1} & {val2:.4f}{str2}\\tabularnewline\n  & ({se1:.4f}) &  ({se2:.4f}) \\tabularnewline\n"""
#table_row_var_2 = """\multirow{{2}}{{*}}{{{var}}} &  & {val2:.4f}{str2} & {val3:.4f}{str3}  \\tabularnewline\n  &  & ({se2:.4f}) & ({se3:.4f}) \\tabularnewline\n"""
#table_row_var_3 = """\multirow{{2}}{{*}}{{{var}}} &  &  & {val3:.4f}{str3}  \\tabularnewline\n  &  &  & ({se3:.4f}) \\tabularnewline\n"""
#table_row_int = '{stat} & {val1:.4g} & {val2:.4g} & {val3:.4g}  \\tabularnewline\n'
#table_row_float = '{stat} & {val1:.3f} & {val2:.3f} & {val3:.3f}  \\tabularnewline\n'
table_row_strings_comp = '{stat} & {val1} & {val2} \\tabularnewline \n'
#table_row_test_1 = '{test} & {val1:.3f} & {val2:.3f} & {val3:.3f}  \\tabularnewline\n'
#table_row_test_2 = '{test} &  & {val2:.3f} & {val3:.3f}  \\tabularnewline\n'
#table_row_test_3 = '{test} &  & & {val3:.3f} \\tabularnewline\n'
table_row_price_comp = """\multirow{{2}}{{*}}{{{var}}} & {val1:.4f}{str1} & {val2:.4f}{str2} \\tabularnewline\n  & ({se1:.4f}) & ({se2:.4f})\\tabularnewline\n"""

# Add significance stars.
eres_base['sig_stars'] = ''
eres_base.loc[eres_base.pval > 0.1, 'sig_stars'] = ''
eres_base.loc[(eres_base.pval <= 0.1) & (
    eres_base.pval > 0.05), 'sig_stars'] = '*'
eres_base.loc[(eres_base.pval <= 0.05) & (
    eres_base.pval > 0.01), 'sig_stars'] = '**'
eres_base.loc[eres_base.pval <= 0.01, 'sig_stars'] = '***'
eres_lsv['sig_stars'] = ''
eres_lsv.loc[eres_lsv.pval > 0.1, 'sig_stars'] = ''
eres_lsv.loc[(eres_lsv.pval <= 0.1) & (
    eres_lsv.pval > 0.05), 'sig_stars'] = '*'
eres_lsv.loc[(eres_lsv.pval <= 0.05) & (
    eres_lsv.pval > 0.01), 'sig_stars'] = '**'
eres_lsv.loc[eres_lsv.pval <= 0.01, 'sig_stars'] = '***'

# Switch order of coefficients and drop meaningless WTPs from column.
# if model_type == 1 or model_type == 7 or model_type == 8:
row_order = [10, 13, 0, 12, 1, 11, 14, 7, 8, 9, 15, 16]
info_order = [7, 8, 9, 15, 16]
# elif model_type == 4:
#     row_order = [11, 14, 0, 13, 1, 12, 15, 16, 7, 8, 9, 10, 17, 18]
#     info_order = [7, 8, 9, 10, 17, 18]
# elif model_type == 5:
#     row_order = [10, 13, 0, 12, 1, 11, 7, 8, 9, 14, 15]
#     info_order = [7, 8, 9, 14, 15]
# elif model_type == 6:
#     row_order = [7, 10, 0, 9, 1, 8, 11, 12, 13]
#     info_order = [12, 13]

# print('Order of original coefficients in table:')
# print(row_order)
# Write the results to a LaTeX table.
with open(project_paths_join('OUT_TABLES', 'estresults_lsv_comp' + '.tex'), 'w') as tex_file:
    # Top of table.
    tex_file.write('\\begin{tabular}{lcc}\n\\toprule\n')
    # Header row with names and dependent variables.
    #     tex_file.write(table_row_strings.format(stat='', val1='Main Model',
    #                                            val2=''))
    tex_file.write(table_row_strings.format(stat='', val1='Baseline',
                                            val2='Larger Logit Smoother'))
    #     tex_file.write(table_row_strings.format(stat='', val1=ols_results['reg_names'][0],
    #                                            val2=ols_results['reg_names'][1],
    #                                            val3=ols_results['reg_names'][2]))
    tex_file.write('\\midrule')
    # legend.at[0, row_idx]
    # Write coefficients to table.
    for i, var in enumerate(row_order):
        # row_idx = row_order[i]
        row_idx = var
        # print(row_idx)
        # if row_idx in(3,4,5,7,9,11,12):
        # if model_type == 1 or model_type == 7 or model_type == 8:
        if row_idx in (10, 13, 11, 7, 8, 9, 15, 16):
            tex_file.write(table_row_price_comp.format(var=eres_base['var_name'][row_idx],
                                                       val1=eres_base['beta'][row_idx],
                                                       str1=eres_base['sig_stars'][row_idx],
                                                       se1=eres_base['sigma'][row_idx],
                                                       val2=eres_lsv['beta'][row_idx],
                                                       str2=eres_lsv['sig_stars'][row_idx],
                                                       se2=eres_lsv['sigma'][row_idx]
                                                       ))
        else:
            tex_file.write(table_row_var_comp.format(var=eres_base['var_name'][row_idx],
                                                     val1=eres_base['beta'][row_idx],
                                                     str1=eres_base['sig_stars'][row_idx],
                                                     se1=eres_base['sigma'][row_idx],
                                                     val2=eres_lsv['beta'][row_idx],
                                                     str2=eres_lsv['sig_stars'][row_idx],
                                                     se2=eres_lsv['sigma'][row_idx]
                                                     ))
    tex_file.write('\\midrule\n')
    # if gmm_type == 1:
    tex_file.write('\\multicolumn{3}{p{11cm}}{\\footnotesize{\\textit{Notes: This table compares the results from our baseline model with a logit smoother value of 0.15 and a robustness check model with a logit smoother value of 0.25. Both models are estimated using GMM with block-diagonal 2SLS weighting matrix. Standard errors in parentheses. *,**,*** denote significance at the 10\%, 5\% and 1\%-level respectively.}}}\\tabularnewline\n')
    # elif gmm_type == 2:
    #     tex_file.write('\\multicolumn{3}{p{11cm}}{\\footnotesize{\\textit{Notes: Results from estimating the demand model using efficient 2-step GMM  weighting matrix. Standard errors in parentheses. *,**,*** denote significance at the 10\%, 5\% and 1\%-level respectively. - denotes non-interpretable willingness-to-pay.}}}\\tabularnewline\n')
    # Bottom of table.
    tex_file.write('\\bottomrule\n\\end{tabular}\n')



####################
# NEW TABLE FOR MAIN TEXT COMBINING BASELINE AND LARGE MODEL RESULTS.

# Update some parts of legend.
legend_4[0][14] = 'Income x price'
legend_4[0][15] = 'Switching cost (non-seniors)'
eres_4['var_name'][14] = 'Income x price'
eres_4['var_name'][15] = 'Switching cost (non-seniors)'

# Set a different types of table row with placeholders.
table_row_strings = '{stat} & {val1} & {val2} & {val3} & {val4} \\tabularnewline \n'

table_row_var_both = """\multirow{{2}}{{*}}{{{var}}} & {val1:.4f}{str1} & {val2: .2f} & {val3:.4f}{str3} & {val4: .2f}  \\tabularnewline\n  & ({se1:.4f}) & & ({se2:.4f}) & \\tabularnewline\n"""
table_row_var_mod_1 = """\multirow{{2}}{{*}}{{{var}}} & {val1:.4f}{str1} & {val2: .2f} &  &  \\tabularnewline\n  & ({se1:.4f}) & &  & \\tabularnewline\n"""
table_row_var_mod_4 = """\multirow{{2}}{{*}}{{{var}}} & &  & {val3:.4f}{str3} & {val4: .2f}  \\tabularnewline\n  & & & ({se2:.4f}) & \\tabularnewline\n"""


table_row_price_both = """\multirow{{2}}{{*}}{{{var}}} & {val1:.4f}{str1} & - & {val2:.4f}{str2} & -  \\tabularnewline\n  & ({se1:.4f}) & & ({se2:.4f}) & \\tabularnewline\n"""
table_row_price_mod_1 = """\multirow{{2}}{{*}}{{{var}}} & {val1:.4f}{str1} & - &  &   \\tabularnewline\n  & ({se1:.4f}) & & & \\tabularnewline\n"""
table_row_price_mod_4 = """\multirow{{2}}{{*}}{{{var}}} &  &  & {val2:.4f}{str2} & -  \\tabularnewline\n  & & & ({se2:.4f}) & \\tabularnewline\n"""


# Switch order of coefficients and drop meaningless WTPs from column.
# if model_type == 1 or model_type == 7 or model_type == 8:
#     row_order = [10, 13, 0, 12, 1, 11, 14, 7, 8, 9, 15, 16]
#     #info_order = [7, 8, 9, 15, 16]
# elif model_type == 4:
row_order = [11, 14, 0, 13, 1, 12, 15, 16, 7, 8, 9, 10, 17, 18]
#info_order = [7, 8, 9, 10, 17, 18]

table_estresults_full_row_header_group = ' & \\multicolumn{{2}}{{c}}{{{val1}}} & \\multicolumn{{2}}{{c}}{{{val2}}} \\tabularnewline \n'

# print('Order of original coefficients in table:')
# print(row_order)
# Write the results to a LaTeX table.
with open(project_paths_join('OUT_TABLES', 'estresults_maintext.tex'), 'w') as tex_file:

    # Top of table.
    tex_file.write('\\begin{tabular}{l|cc|cc}\n\\toprule\n')
    # Header row with names and dependent variables.
    tex_file.write(table_estresults_full_row_header_group.format(val1=r'\textbf{Baseline}',val2=r'\textbf{Extension}'))
    tex_file.write(table_row_strings.format(stat='', val1='Coefficients',
                                            val2='WTP in EUR', val3='Coefficients',
                                            val4='WTP in EUR'))

    tex_file.write('\\midrule')
# legend.at[0, row_idx]
    # Write coefficients to table.
    for i, var in enumerate(row_order):
        # row_idx = row_order[i]
        row_idx = var
        # print(row_idx)
        # if row_idx in(3,4,5,7,9,11,12):
        # Adjust index for model 1 spec.
        if row_idx<10:
            mod_1_idx = row_idx
        elif row_idx>=10 and row_idx<17:
            mod_1_idx = row_idx - 1
        elif row_idx>=17:
            mod_1_idx = row_idx - 2

        if row_idx in (11, 14, 12, 7, 8, 9, 10, 17, 18):
            if row_idx==10:
                tex_file.write(table_row_price_mod_4.format(var=eres_4['var_name'][row_idx], val2=eres_4['beta'][row_idx],
                                                    str2=eres_4['sig_stars'][row_idx],
                                                    se2=eres_4['sigma'][row_idx]))    
            else:
                tex_file.write(table_row_price_both.format(var=eres_4['var_name'][row_idx], val1=eres_1['beta'][mod_1_idx],
                                                        str1=eres_1['sig_stars'][mod_1_idx],
                                                        se1=eres_1['sigma'][mod_1_idx],
                                                        val2=eres_4['beta'][row_idx],
                                                        str2=eres_4['sig_stars'][row_idx],
                                                        se2=eres_4['sigma'][row_idx]
                                                        ))
        else:
            if row_idx==15:
                # Add additional line for switching cost estimates in homogeneous   baseline spec.
                tex_file.write(table_row_var_mod_1.format(var=eres_1['var_name'][mod_1_idx], val1=eres_1['beta'][mod_1_idx],
                                                val2=eres_1['WTP'][mod_1_idx],
                                                str1=eres_1['sig_stars'][mod_1_idx],
                                                se1=eres_1['sigma'][mod_1_idx]))            
                # Add additional line for non-senior switching costs in baseline spec.
                tex_file.write(table_row_var_mod_4.format(var=eres_4['var_name'][row_idx], val3=eres_4['beta'][row_idx],
                                                val4=eres_4['WTP'][row_idx],
                                                str3=eres_4['sig_stars'][row_idx],
                                                se2=eres_4['sigma'][row_idx]))            

            elif row_idx==16:
                # Add additional line for senior switching cost in extended spec.
                tex_file.write(table_row_var_mod_4.format(var=eres_4['var_name'][row_idx], val3=eres_4['beta'][row_idx],
                                                val4=eres_4['WTP'][row_idx],
                                                str3=eres_4['sig_stars'][row_idx],
                                                se2=eres_4['sigma'][row_idx]))  
            else:
                tex_file.write(table_row_var_both.format(var=eres_4['var_name'][row_idx], val1=eres_1['beta'][mod_1_idx],
                                                    val2=eres_1['WTP'][mod_1_idx],
                                                    str1=eres_1['sig_stars'][mod_1_idx],
                                                    se1=eres_1['sigma'][mod_1_idx],
                                                    val3=eres_4['beta'][row_idx],
                                                    val4=eres_4['WTP'][row_idx],
                                                    str3=eres_4['sig_stars'][row_idx],
                                                    se2=eres_4['sigma'][row_idx]))
    tex_file.write('\\midrule\n')
    tex_file.write('\\multicolumn{5}{p{16cm}}{\\footnotesize{\\textit{Notes: Results from estimating the demand model using efficient 2-step GMM  weighting matrix. The left panel contains our baseline specification with homogeneous PCW search costs and switching costs. The right panel contains our model extension in which PCW search costs and switching costs vary across age groups (non-seniors and seniors). Both specifications include firm fixed effects. Standard errors in parentheses. *,**,*** denote significance at the 10\%, 5\% and 1\%-level respectively. - denotes non-interpretable willingness-to-pay. }}}\\tabularnewline\n')
    # Bottom of table.
    tex_file.write('\\bottomrule\n\\end{tabular}\n')

################################################################################
# Main estimation results table in AEJ:Micro format.
table_estresults_full_row_header_group = ' & \\multicolumn{{2}}{{c}}{{{val1}}} & \\multicolumn{{2}}{{c}}{{{val2}}} \\tabularnewline \n'
# Important: Don't display significance stars!
# Write the results to a LaTeX table.
with open(project_paths_join('OUT_TABLES', 'estresults_maintext_aej.tex'), 'w') as tex_file:

    # Top of table.
    tex_file.write('\\begin{tabular}{lcccc}\n\\toprule\n')
    # Header row with names and dependent variables.
    tex_file.write(table_estresults_full_row_header_group.format(val1=r'\textbf{Baseline}',val2=r'\textbf{Extension}'))
    tex_file.write(table_row_strings.format(stat='', val1='Coefficients',
                                            val2='WTP in EUR', val3='Coefficients',
                                            val4='WTP in EUR'))

    tex_file.write('\\midrule')
# legend.at[0, row_idx]
    # Write coefficients to table.
    for i, var in enumerate(row_order):
        # row_idx = row_order[i]
        row_idx = var
        # print(row_idx)
        # if row_idx in(3,4,5,7,9,11,12):
        # Adjust index for model 1 spec.
        if row_idx<10:
            mod_1_idx = row_idx
        elif row_idx>=10 and row_idx<17:
            mod_1_idx = row_idx - 1
        elif row_idx>=17:
            mod_1_idx = row_idx - 2

        if row_idx in (11, 14, 12, 7, 8, 9, 10, 17, 18):
            if row_idx==10:
                tex_file.write(table_row_price_mod_4.format(var=eres_4['var_name'][row_idx], val2=eres_4['beta'][row_idx],
                                                    str2='',
                                                    se2=eres_4['sigma'][row_idx]))    
            else:
                tex_file.write(table_row_price_both.format(var=eres_4['var_name'][row_idx], val1=eres_1['beta'][mod_1_idx],
                                                        str1='',
                                                        se1=eres_1['sigma'][mod_1_idx],
                                                        val2=eres_4['beta'][row_idx],
                                                        str2='',
                                                        se2=eres_4['sigma'][row_idx]
                                                        ))
        else:
            if row_idx==15:
                # Add additional line for switching cost estimates in homogeneous   baseline spec.
                tex_file.write(table_row_var_mod_1.format(var=eres_1['var_name'][mod_1_idx], val1=eres_1['beta'][mod_1_idx],
                                                val2=eres_1['WTP'][mod_1_idx],
                                                str1='',
                                                se1=eres_1['sigma'][mod_1_idx]))            
                # Add additional line for non-senior switching costs in baseline spec.
                tex_file.write(table_row_var_mod_4.format(var=eres_4['var_name'][row_idx], val3=eres_4['beta'][row_idx],
                                                val4=eres_4['WTP'][row_idx],
                                                str3='',
                                                se2=eres_4['sigma'][row_idx]))            

            elif row_idx==16:
                # Add additional line for senior switching cost in extended spec.
                tex_file.write(table_row_var_mod_4.format(var=eres_4['var_name'][row_idx], val3=eres_4['beta'][row_idx],
                                                val4=eres_4['WTP'][row_idx],
                                                str3='',
                                                se2=eres_4['sigma'][row_idx]))  
            else:
                tex_file.write(table_row_var_both.format(var=eres_4['var_name'][row_idx], val1=eres_1['beta'][mod_1_idx],
                                                    val2=eres_1['WTP'][mod_1_idx],
                                                    str1='',
                                                    se1=eres_1['sigma'][mod_1_idx],
                                                    val3=eres_4['beta'][row_idx],
                                                    val4=eres_4['WTP'][row_idx],
                                                    str3='',
                                                    se2=eres_4['sigma'][row_idx]))
    # tex_file.write('\\midrule\n')
    # tex_file.write('\\multicolumn{5}{p{16cm}}{\\footnotesize{\\textit{Notes: Results from estimating the demand model using efficient 2-step GMM  weighting matrix. The left panel contains our baseline specification with homogeneous PCW search costs and switching costs. The right panel contains our model extension in which PCW search costs and switching costs vary across age groups (non-seniors and seniors). Both specifications include firm fixed effects. Standard errors in parentheses. *,**,*** denote significance at the 10\%, 5\% and 1\%-level respectively. - denotes non-interpretable willingness-to-pay. }}}\\tabularnewline\n')
    # Bottom of table.
    tex_file.write('\\bottomrule\n\\end{tabular}\n')


################################################################################
# NEW TABLE FOR REF DOC TO COMPARE BASELINE WITH AND WITHOUT HAUSMAN IV.

# Update some parts of legend.
legend_4[0][13] = 'Income x price'
# legend_4[0][14] = 'Switching cost'
eres_1['var_name'][13] = 'Income x price'
# eres_1['var_name'][14] = 'Switching cost'

# Set a different types of table row with placeholders.
table_row_strings = '{stat} & {val1} & {val2} & {val3} & {val4} \\tabularnewline \n'

table_row_var_both = """\multirow{{2}}{{*}}{{{var}}} & {val1:.4f}{str1} & {val2: .2f} & {val3:.4f}{str3} & {val4: .2f}  \\tabularnewline\n  & ({se1:.4f}) & & ({se2:.4f}) & \\tabularnewline\n"""
# table_row_var_mod_1 = """\multirow{{2}}{{*}}{{{var}}} & {val1:.4f}{str1} & {val2: .2f} &  &  \\tabularnewline\n  & ({se1:.4f}) & &  & \\tabularnewline\n"""
# table_row_var_mod_4 = """\multirow{{2}}{{*}}{{{var}}} & &  & {val3:.4f}{str3} & {val4: .2f}  \\tabularnewline\n  & & & ({se2:.4f}) & \\tabularnewline\n"""
table_row_price_both = """\multirow{{2}}{{*}}{{{var}}} & {val1:.4f}{str1} & - & {val2:.4f}{str2} & -  \\tabularnewline\n  & ({se1:.4f}) & & ({se2:.4f}) & \\tabularnewline\n"""
# table_row_price_mod_1 = """\multirow{{2}}{{*}}{{{var}}} & {val1:.4f}{str1} & - &  &   \\tabularnewline\n  & ({se1:.4f}) & & & \\tabularnewline\n"""
# table_row_price_mod_4 = """\multirow{{2}}{{*}}{{{var}}} &  &  & {val2:.4f}{str2} & -  \\tabularnewline\n  & & & ({se2:.4f}) & \\tabularnewline\n"""


# Switch order of coefficients and drop meaningless WTPs from column.
# Model 1 and model 7 have same ordering.
row_order = [10, 13, 0, 12, 1, 11, 14, 7, 8, 9, 15, 16]

# print('Order of original coefficients in table:')
# print(row_order)
# Write the results to a LaTeX table.
with open(project_paths_join('OUT_TABLES', 'estresults_hausmancomp.tex'), 'w') as tex_file:

    # Top of table.
    tex_file.write('\\begin{tabular}{l|cc|cc}\n\\toprule\n')
    # Header row with names and dependent variables.
    tex_file.write(table_row_strings.format(stat='', val1=r'\textbf{Baseline}',
                                            val2='',val3=r'\textbf{No Hausman IV}',val4=''))
    tex_file.write(table_row_strings.format(stat='', val1='Coefficients',
                                            val2='WTP in EUR', val3='Coefficients',
                                            val4='WTP in EUR'))

    tex_file.write('\\midrule')
# legend.at[0, row_idx]
    # Write coefficients to table.
    for i, var in enumerate(row_order):
        # row_idx = row_order[i]
        row_idx = var
        mod_1_idx = row_idx
        # print(row_idx)
        # if row_idx in(3,4,5,7,9,11,12):
        # Adjust index for model 1 spec.
        # if row_idx<10:
        #     mod_1_idx = row_idx
        # elif row_idx>=10 and row_idx<17:
        #     mod_1_idx = row_idx - 1
        # elif row_idx>=17:
        #     mod_1_idx = row_idx - 2

        if row_idx in (10, 13, 11, 7, 8, 9, 15, 16):
            tex_file.write(table_row_price_both.format(var=eres_7['var_name'][row_idx], val1=eres_1['beta'][mod_1_idx],
                                                        str1=eres_1['sig_stars'][mod_1_idx],
                                                        se1=eres_1['sigma'][mod_1_idx],
                                                        val2=eres_7['beta'][row_idx],
                                                        str2=eres_7['sig_stars'][row_idx],
                                                        se2=eres_7['sigma'][row_idx]
                                                        ))
        else:
            tex_file.write(table_row_var_both.format(var=eres_7['var_name'][row_idx], val1=eres_1['beta'][mod_1_idx],
                                                    val2=eres_1['WTP'][mod_1_idx],
                                                    str1=eres_1['sig_stars'][mod_1_idx],
                                                    se1=eres_1['sigma'][mod_1_idx],
                                                    val3=eres_7['beta'][row_idx],
                                                    val4=eres_7['WTP'][row_idx],
                                                    str3=eres_7['sig_stars'][row_idx],
                                                    se2=eres_7['sigma'][row_idx]))
    tex_file.write('\\midrule\n')
    tex_file.write('\\multicolumn{5}{p{16cm}}{\\footnotesize{\\textit{Notes: Results from estimating the demand model using efficient 2-step GMM  weighting matrix. The left panel contains our baseline specification that uses our Hausman IV. The right panel contains the same model specification, but does not use the Hausman IV. Both specifications include firm fixed effects. Standard errors in parentheses. *,**,*** denote significance at the 10\%, 5\% and 1\%-level respectively. - denotes non-interpretable willingness-to-pay. }}}\\tabularnewline\n')
    # Bottom of table.
    tex_file.write('\\bottomrule\n\\end{tabular}\n')