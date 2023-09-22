"""
    Write several tables and graphs with summary of counterfactual results into tex-table.
    Plot results for all of the different model types.

    1. Plot counterfactual evolution of PCW usage.
    2. Format comparison of market shares table.
    3. Counterfactual results table comparing model 1 and 4 with everything.
    4. Compare selected welfare statistics for full and restricted models.
    5. Plot indifference-markup price paths for different specifications.

    OUTPUT TABLES (one separate file for each combination of model_type and GMM type, unless otherwise indicated):
    - cf_full_gmm_GMMTYPE.tex (comparison across the two main models)
    - cf_comp_gmm_GMMTYPE.tex (comparison across the four main models)
    - cf_shares_FILESUFFIX.tex 
    - cf_shares_small_FILESUFFIX.tex (probably not used in paper)

    OUTPUT FIGURES (one separate file for each combination of model_type and GMM type, unless otherwise indicated):
    - pcw_usage_obspred_FILESUFFIX.pdf (not used in paper)
    - cf_pcw_usage_FILESUFFIX.pdf (not used in paper)
    - cf_pcw_usage_comb_FILESUFFIX.pdf
    - cf_pcw_change_FILESUFFIX.pdf (not used in paper)
    - cf_price_id_FILESUFFIX.pdf

    Most importantly, this file formats Table 3 in the main text.

"""
from bld.project_paths import project_paths_join
import numpy as np
import pandas as pd
#import scipy.io as sio
import matplotlib.dates as mdates
from matplotlib import pyplot
pd.options.mode.chained_assignment = None  # default='warn'


model_type = 2
gmm_type = 2
gmm_type_restrict = 1
file_suffix = str(gmm_type) + '_mod' + str(model_type)
file_suffix_modelonly = '1' + '_mod' + str(model_type)

# Create overview table for all counterfactual results using both full models.

# Load data from both model 1 and model 4.
cs_full_mod1 = pd.read_csv(project_paths_join(
    'OUT_ANALYSIS', 'welfare_comp_' + str(gmm_type) + '_mod1.csv'))
# Adjust label for regulator campaign to reflect correct units.
cs_full_mod1['Row'][2] = 'Info campaign (total gain in mio. EUR)'
cs_full_mod1['Row'][3] = 'Automatic lowest price'
cs_full_mod1 = cs_full_mod1.set_index('Row')

cs_full_mod4 = pd.read_csv(project_paths_join(
    'OUT_ANALYSIS', 'welfare_comp_' + str(gmm_type) + '_mod4.csv'))
cs_full_mod4['Row'][2] = 'Info campaign (total gain in mio. EUR)'
cs_full_mod4['Row'][3] = 'Automatic lowest price'
cs_full_mod4 = cs_full_mod4.set_index('Row')
# Load data for restricted models.
cs_full_mod5 = pd.read_csv(project_paths_join(
    'OUT_ANALYSIS', 'welfare_comp_' + str(gmm_type_restrict) + '_mod5.csv'))
cs_full_mod5['Row'][2] = 'Info campaign (total gain in mio. EUR)'
cs_full_mod5['Row'][3] = 'Automatic lowest price'
cs_full_mod5 = cs_full_mod5.set_index('Row')

cs_full_mod6 = pd.read_csv(project_paths_join(
    'OUT_ANALYSIS', 'welfare_comp_' + str(gmm_type_restrict) + '_mod6.csv'))
cs_full_mod6['Row'][2] = 'Info campaign (total gain in mio. EUR)'
cs_full_mod6['Row'][3] = 'Automatic lowest price'
cs_full_mod6 = cs_full_mod6.set_index('Row')


# Share of seniors for weighting welfare gains from regulator campaign.
share_senior = 0.19
# Correct values for total welfare gain from regulator campaign of different sub groups.
cs_full_mod1['Non-seniors'][2] = cs_full_mod1['Non-seniors'][2] * \
    (1-share_senior)
cs_full_mod1['Seniors'][2] = cs_full_mod1['Seniors'][2] * share_senior
cs_full_mod4['Non-seniors'][2] = cs_full_mod4['Non-seniors'][2] * \
    (1-share_senior)
cs_full_mod4['Seniors'][2] = cs_full_mod4['Seniors'][2] * share_senior
cs_full_mod5['Non-seniors'][2] = cs_full_mod5['Non-seniors'][2] * \
    (1-share_senior)
cs_full_mod5['Seniors'][2] = cs_full_mod5['Seniors'][2] * share_senior
cs_full_mod6['Non-seniors'][2] = cs_full_mod6['Non-seniors'][2] * \
    (1-share_senior)
cs_full_mod6['Seniors'][2] = cs_full_mod6['Seniors'][2] * share_senior
# Double check that calculations are correct.
test_rc_calc = cs_full_mod1['Non-seniors'][2] + cs_full_mod1['Seniors'][2]


# Setup table templates for full counterfactual table with detailed results.
table_cf_full_row_header_group = ' & \\multicolumn{{2}}{{c}}{{{val1}}} & \\multicolumn{{2}}{{c}}{{{val2}}} &  \\multicolumn{{2}}{{c}}{{{val3}}} \\tabularnewline \n'
table_cf_full_row_header_mod = ' (in EUR per month and consumer) & {val1} & {val2} & {val1} & {val2} & {val1} & {val2}  \\tabularnewline \n'
table_cf_full_row = '{cf_name} & {val1:.2f} & {val2:.2f} & {val3:.2f} & {val4:.2f} & {val5: .2f} & {val6: .2f} \\tabularnewline \n'

# Write to LaTeX table: Large table with observed vs. predicted
with open(project_paths_join('OUT_TABLES', 'cf_full_gmm_' + str(gmm_type)+'.tex'), 'w') as tex_file:
    # Top of table.
    tex_file.write('\\begin{tabular}{r|cc|cc|cc}\n\\toprule\n')
    # Header row with names and dependent variables.
    tex_file.write(table_cf_full_row_header_group.format(val1='Average/Total',
                                                         val2='Non-seniors',
                                                         val3='Seniors'))
    tex_file.write(table_cf_full_row_header_mod.format(val1='Baseline',
                                                       val2='Extension'))
    tex_file.write('\\midrule \n')
    # Write counterfactual results to table by iterating over rows from model 1.
    for index, row in cs_full_mod1.iterrows():
        if index!='Incumbent only conv.':
            tex_file.write(table_cf_full_row.format(cf_name=index,
                                                    val1=cs_full_mod1['Average'][index],
                                                    val2=cs_full_mod4['Average'][index],
                                                    val3=cs_full_mod1['Non-seniors'][index],
                                                    val4=cs_full_mod4['Non-seniors'][index],
                                                    val5=cs_full_mod1['Seniors'][index],
                                                    val6=cs_full_mod4['Seniors'][index]))
        if index in ['Reduced PCW search costs', 'Automatic lowest price']:
            tex_file.write('\\midrule  \n')
    # Write bottom of table and notes.
    tex_file.write('\\midrule  \n')
    tex_file.write('\\multicolumn{7}{p{17.5cm}}{\\footnotesize{\\textit{Notes: The table summarizes the welfare results from all couterfactual simulations for both the baseline model with homogeneous market frictions and the model extension with heterogeneous PCW search cost and switching cost. For the counterfactual "Info Campaign" the number indicates the total welfare gain for the whole sample period and the whole region of Flanders. For all other counterfactuals, the numbers indicate the welfare gain in EUR/month for the average consumer. Row \emph{Incumbent all (no SC)} indicates the best-case scenario in which one firm offers all currently available contracts without any market frictions.}}}\\tabularnewline\n')
    # Bottom of table.
    tex_file.write('\\bottomrule\n\\end{tabular}\n')

################################################################################
# Counterfactual table with formatting for AEJ:Micro
with open(project_paths_join('OUT_TABLES', 'cf_full_gmm_' + str(gmm_type)+'_aej'+'.tex'), 'w') as tex_file:
    # Top of table.
    tex_file.write('\\begin{tabular}{rcccccc}\n\\toprule\n')
    # Header row with names and dependent variables.
    tex_file.write(table_cf_full_row_header_group.format(val1='Average/Total',
                                                         val2='Non-seniors',
                                                         val3='Seniors'))
    tex_file.write(table_cf_full_row_header_mod.format(val1='Baseline',
                                                       val2='Extension'))
    tex_file.write('\\midrule \n')
    # Write counterfactual results to table by iterating over rows from model 1.
    for index, row in cs_full_mod1.iterrows():
        if index!='Incumbent only conv.':
            tex_file.write(table_cf_full_row.format(cf_name=index,
                                                    val1=cs_full_mod1['Average'][index],
                                                    val2=cs_full_mod4['Average'][index],
                                                    val3=cs_full_mod1['Non-seniors'][index],
                                                    val4=cs_full_mod4['Non-seniors'][index],
                                                    val5=cs_full_mod1['Seniors'][index],
                                                    val6=cs_full_mod4['Seniors'][index]))
        if index in ['Reduced PCW search costs', 'Automatic lowest price']:
            tex_file.write('\\midrule  \n')
    # Write bottom of table and notes.
    # tex_file.write('\\midrule  \n')
    # tex_file.write('\\multicolumn{7}{p{19.5cm}}{\\footnotesize{\\textit{Notes: The table summarizes the welfare results from all couterfactual simulations for both the baseline model with homogeneous market frictions and the large model with heterogeneous PCW search cost and switching cost. For the counterfactual "Info Campaign" the number indicates the total welfare gain for the whole sample period and the whole region of Flanders. For all other counterfactuals, the numbers indicate the welfare gain in EUR/month for the average consumer. Row \emph{Incumbent all (no SC)} indicates the best-case scenario in which one firm offers all currently available contracts without any market frictions.}}}\\tabularnewline\n')
    # Bottom of table.
    tex_file.write('\\bottomrule\n\\end{tabular}\n')

# Create a smaller comparison table that compares most important counterfactual results for the 2 full and the two restricted models.
# Setup table templates for full counterfactual table with detailed results.
# table_cf_comp_row_header_group = ' & {val1} & {val2} & &  {val3} & \\tabularnewline \n'
table_cf_comp_row_header_mod = ' & {val1} & {val2} & {val3} & {val4} \\tabularnewline \n'
table_cf_comp_row = '{cf_name} & {val1:.2f} & {val2:.2f} & {val3:.2f} & {val4:.2f} \\tabularnewline \n'
cf_comp_list = ['Reduced switching costs',
                'Reduced PCW search costs', 'Info campaign', 'Automatic lowest price',
                'Incumbent only', 'Incumbent only (no SC)', 
                'Incumbent only (logit corr.)', 'Incumbent only (logit corr. and no SC)', 'Incumbent all (no SC)']
# Write to LaTeX table: Large table with observed vs. predicted
with open(project_paths_join('OUT_TABLES', 'cf_comp_gmm_' + str(gmm_type)+'.tex'), 'w') as tex_file:
    # Top of table.
    tex_file.write('\\begin{tabular}{r|cccc}\n\\toprule\n')
    # Header row with names and dependent variables.
    tex_file.write(table_cf_comp_row_header_mod.format(val1='Baseline',
                                                       val2='Extension',
                                                       val3=r'No switching costs',
                                                       val4=r'No PCW search costs'))
    tex_file.write('\\midrule \n')
    # Write counterfactual results to table by iterating over rows from model 1.
    for index, row in cs_full_mod1.iterrows():
        # print(index)
        if index in cf_comp_list:

            if index in ['Reduced switching costs']:
                tex_file.write(table_cf_comp_row.format(cf_name=index,
                                                        val1=cs_full_mod1['Average'][index],
                                                        val2=cs_full_mod4['Average'][index],
                                                        val3=0.0,
                                                        val4=cs_full_mod6['Average'][index]))
            elif index in ['Reduced PCW search costs']:
                tex_file.write(table_cf_comp_row.format(cf_name=index,
                                                        val1=cs_full_mod1['Average'][index],
                                                        val2=cs_full_mod4['Average'][index],
                                                        val3=cs_full_mod5['Average'][index],
                                                        val4=0.0))
            else:
                tex_file.write(table_cf_comp_row.format(cf_name=index,
                                                        val1=cs_full_mod1['Average'][index],
                                                        val2=cs_full_mod4['Average'][index],
                                                        val3=cs_full_mod5['Average'][index],
                                                        val4=cs_full_mod6['Average'][index]))

        if index in ['Reduced PCW search costs', 'Automatic lowest price']:
            tex_file.write('\\midrule  \n')
    # Write bottom of table and notes.
    tex_file.write('\\midrule  \n')
    tex_file.write(r'\multicolumn{5}{p{18.5cm}}{\footnotesize{\textit{Notes: The table compares the welfare results from several couterfactual simulations for both the baseline model with homogeneous market frictions and the model extension with heterogeneous PCW search cost and switching cost, as well as the restricted models that either ignore switching costs ($\psi$) or PCW search costs ($\kappa$). For all counterfactuals, the numbers indicate the welfare gain in EUR/month for the average consumer. Row \emph{Incumbent all (no SC)} indicates the best-case scenario in which one firm offers all currently available contracts without any market frictions.}}}\tabularnewline')
    #Combinations of model type and counterfactual that do not make sense are indicated with $0.0$.
    # Bottom of table.
    tex_file.write('\n\\bottomrule\n\\end{tabular}\n')


# Select model types & GMM type to format.
model_list = [1, 4, 5, 6]
gmm_list = [
    [2],
    [2],
    [1],
    [1]]
model_print = r'Formatting estimation results for model {mt:d} and GMM stage {gmm:d}...\n'
markup_print = r"""Indifference markups for model {mt:d} and GMM stage {gmm:d} (in percent):
Inc only conv: {ioc:.2f}
Inc only: {io:.2f}
Inc only (no SC): {ionosc:.2f}
Inc only (LC, no SC): {iolcnosc:.2f}
Inc all (no SC): {incallnosc:.2f}
Observed: {obs:.2f}

"""

# Initiate empty data frame for filling with contract valuations for comparison across model types.
model_names = ["Baseline", "Large", "No $\psi$", "No $\kappa$"]
# Data frame for comparing different contract valuations implied by different model specifications. Two data frames for two different GMM stages.
comp_contract_val_1 = pd.DataFrame(0, index=np.arange(11), columns=model_names)
comp_contract_val_2 = pd.DataFrame(0, index=np.arange(11), columns=model_names)
# For model comparison, we print 7 estimation statistics
comp_est_stats_1 = pd.DataFrame(0, index=np.arange(7), columns=model_names)
comp_est_stats_2 = pd.DataFrame(0, index=np.arange(7), columns=model_names)


for iter, model_type in enumerate(model_list):
    #print(iter)
    for gmm_type in gmm_list[iter]:
        print(model_print.format(mt=model_type, gmm=gmm_type))
        file_suffix = str(gmm_type) + '_mod' + str(model_type)
        file_suffix_modelonly = '1' + '_mod' + str(model_type)

        # Load indifference markups.
        id_mu_data = pd.read_csv(project_paths_join(
            'OUT_ANALYSIS', 'id_mu_' + file_suffix + '.csv'))
        id_mu_data = id_mu_data.set_index('Row')

        pcw_comp = pd.read_csv(project_paths_join(
            'OUT_ANALYSIS', 'cf_pcw_use_' + file_suffix+'.csv'))
        pcw_diff = pd.read_csv(project_paths_join(
            'OUT_ANALYSIS', 'cf_pcw_change_'+file_suffix+'.csv'))
        # Label columns.
        pcw_comp.columns = ['Observed', 'Predicted',
                            'Switching Subsidy', 'Search Cost Reduction', 'Info campagin']
        pcw_diff.columns = ['Switching Subsidy',
                            'Search Cost Reduction', 'Info campaign']
        # Construct data axis.
        date_num = mdates.datestr2num(['2012-02', '2012-03', '2012-04', '2012-05', '2012-06',
                                       '2012-07', '2012-08', '2012-09', '2012-10', '2012-11', '2012-12',
                                       '2013-01', '2013-02', '2013-03', '2013-04', '2013-05', '2013-06',
                                       '2013-07', '2013-08', '2013-09', '2013-10', '2013-11', '2013-12',
                                       '2014-01', '2014-02', '2014-03', '2014-04', '2014-05', '2014-06',
                                       '2014-07', '2014-08', '2014-09', '2014-10', '2014-11', '2014-12',
                                       '2015-01', '2015-02', '2015-03', '2015-04', '2015-05', '2015-06',
                                       '2015-07', '2015-08', '2015-09', '2015-10', '2015-11', '2015-12',
                                       '2016-01', '2016-02', '2016-03', '2016-04', '2016-05', '2016-06'])

        # Define the basic characteristics of the figure.
        # Observed vs. predicted PCW usage.
        pyplot.clf()
        fig, ax = pyplot.subplots(figsize=(8, 6))
        ax = pyplot.subplot(111)
        pyplot.grid(True)
        pyplot.title('Aggregate PCW usage over time')
        # Plot price index datafor each operator.
        pyplot.plot_date(date_num, pcw_comp['Observed'],
                         'k--', label='Data', xdate=True, ydate=False)
        pyplot.plot_date(date_num, pcw_comp['Predicted'],
                         'r--', label='Predicted',  xdate=True, ydate=False)
        ax.set_ylabel('Monthly share of PCW users')
        ax.set_ylim([0, 0.2])
        pyplot.xticks(np.arange(min(date_num), max(date_num)+1, 350))
        # Shink current axis by 17%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + 0.13, box.width, box.height * 0.83])
        ax.set_xlim([mdates.datestr2num('2012-02'),
                     mdates.datestr2num('2016-06')])
        # pyplot.xticks(np.arange(min(date_num), max(date_num)+180, 365.25))
        ax.xaxis.set_major_locator(mdates.MonthLocator(
            interval=10))  # to get a tick every 15 minutes
        ax.xaxis.set_major_formatter(
            mdates.DateFormatter('%b-%Y'))  # optional formatting
        # Put a legend to the right of the current axis
        pyplot.legend(loc='upper center', ncol=2)
        ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.28), ncol=2)
        # # Save the figure.
        pyplot.savefig(project_paths_join(
            'OUT_FIGURES', 'pcw_usage_obspred_'+file_suffix+'.pdf'), bbox_inches='tight')

        # Counterfactual PCW usage.
        pyplot.clf()
        fig, ax = pyplot.subplots(figsize=(8, 6))
        ax = pyplot.subplot(111)
        pyplot.grid(True)
        pyplot.title('Counterfactual PCW usage over time')
        pyplot.plot_date(date_num, pcw_comp['Switching Subsidy'], 'g-',
                         label='Reduced switching cost', xdate=True, ydate=False)
        pyplot.plot_date(date_num, pcw_comp['Search Cost Reduction'],
                         'y-', label='Reduced PCW search cost', xdate=True, ydate=False)
        ax.set_ylabel('Monthly share of PCW users')
        ax.set_ylim([0, 0.7])
        pyplot.xticks(np.arange(min(date_num), max(date_num)+1, 350))
        # Shink current axis by 17%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + 0.13, box.width, box.height * 0.83])
        ax.set_xlim([mdates.datestr2num('2012-02'),
                     mdates.datestr2num('2016-06')])
        # pyplot.xticks(np.arange(min(date_num), max(date_num)+180, 365.25))
        ax.xaxis.set_major_locator(mdates.MonthLocator(
            interval=10))  # to get a tick every 15 minutes
        ax.xaxis.set_major_formatter(
            mdates.DateFormatter('%b-%Y'))  # optional formatting
        # Put a legend to the right of the current axis
        pyplot.legend(loc='upper center', ncol=2)
        ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.28), ncol=2)
        # # Save the figure.
        pyplot.savefig(project_paths_join(
            'OUT_FIGURES', 'cf_pcw_usage_' + file_suffix+'.pdf'), bbox_inches='tight')

        # Compare observed and counterfactual PCW usage.
        # Counterfactual PCW usage.
        pyplot.clf()
        fig, ax = pyplot.subplots(figsize=(8, 6))
        ax = pyplot.subplot(111)
        pyplot.grid(True)
        #pyplot.title('Observed and counterfactual PCW usage over time')
        # Plot price index datafor each operator.
        pyplot.plot_date(date_num, pcw_comp['Observed'],
                         'k--', label='Data', xdate=True, ydate=False)
        pyplot.plot_date(date_num, pcw_comp['Predicted'],
                         'r--', label='Predicted',  xdate=True, ydate=False)
        pyplot.plot_date(date_num, pcw_comp['Switching Subsidy'], 'g-',
                         label='Reduced switching cost', xdate=True, ydate=False)
        pyplot.plot_date(date_num, pcw_comp['Search Cost Reduction'],
                         'y-', label='Reduced PCW search cost', xdate=True, ydate=False)
        ax.set_ylabel('Monthly share of PCW users')
        ax.set_ylim([0, 0.8])
        pyplot.xticks(np.arange(min(date_num), max(date_num)+1, 350))
        # Shink current axis by 17%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + 0.13, box.width, box.height * 0.83])
        ax.set_xlim([mdates.datestr2num('2012-02'),
                     mdates.datestr2num('2016-06')])
        # pyplot.xticks(np.arange(min(date_num), max(date_num)+180, 365.25))
        ax.xaxis.set_major_locator(mdates.MonthLocator(
            interval=10))  # to get a tick every 15 minutes
        ax.xaxis.set_major_formatter(
            mdates.DateFormatter('%b-%Y'))  # optional formatting
        # Put a legend to the right of the current axis
        pyplot.legend(loc='upper center', ncol=2)
        ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.28), ncol=3)
        # # Save the figure.
        pyplot.savefig(project_paths_join(
            'OUT_FIGURES', 'cf_pcw_usage_comb_' + file_suffix + '.pdf'), bbox_inches='tight')

        # Same for changes in PCW usage.
        # Define the basic characteristics of the figure.
        pyplot.clf()
        fig, ax = pyplot.subplots(figsize=(8, 6))
        # pyplot.plot()
        pyplot.grid(True)
        # pyplot.ylabel('Market shares in %')
        #pyplot.title('Counterfactuals: Changes in Aggregate PCW Usage')
        # ax = pyplot.subplot(111)
        # Plot price index datafor each operator.
        # pyplot.plot_date(date_num, pcw_comp['Observed'] , 'k--',label='Data', xdate=True,ydate=False)
        # pyplot.plot_date(date_num, pcw_comp['Predicted'],'r-',label='Predicted',  xdate=True,ydate=False)
        pyplot.plot_date(date_num, pcw_diff['Switching Subsidy'], 'r-',
                         label='Reduced switching cost', xdate=True, ydate=False)
        pyplot.plot_date(date_num, pcw_diff['Search Cost Reduction'],
                         'b-', label='Reduced PCW search cost', xdate=True, ydate=False)
        # pyplot.plot_date(date_num, ms_obs_lampiris,'c-',label='Lampiris', xdate=True,ydate=False)
        # pyplot.plot_date(date_num, ms_obs_other, 'k--',label='Other', xdate=True,ydate=False)
        ax.set_ylabel('Change in share of PCW users (in %-points)')
        ax.set_ylim([0.1, 0.6])
        ax.set_xlim([mdates.datestr2num('2012-02'),
                     mdates.datestr2num('2016-06')])
        # pyplot.xticks(np.arange(min(date_num), max(date_num)+180, 365.25))
        ax.xaxis.set_major_locator(mdates.MonthLocator(
            interval=10))  # to get a tick every 15 minutes
        ax.xaxis.set_major_formatter(
            mdates.DateFormatter('%b-%Y'))  # optional formatting
        # Shink current axis by 17%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + 0.13, box.width, box.height * 0.83])
        ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=2)
        # Save the figure.
        pyplot.savefig(project_paths_join(
            'OUT_FIGURES', 'cf_pcw_change_' + file_suffix + '.pdf'), bbox_inches='tight')

        ################################################################################
        # Format market share table for counterfactual summary.
        # Load estimation results for gross and net sample from csv-file.
        shares = pd.read_csv(project_paths_join(
            'OUT_ANALYSIS', 'cf_marketshares_'+file_suffix+'.csv'))

        # # Load row and column labels.
        # row_labels = pd.read_csv(project_paths_join(
        #     'OUT_ANALYSIS', 'cfms_rowlabels_' + file_suffix+'.csv'), header=None)
        # col_labels = pd.read_csv(project_paths_join(
        #     'OUT_ANALYSIS', 'cfms_collabels_' + file_suffix+'.csv'))

        # Rename variables nicely.
        shares.columns = ['Observed', 'Predicted',
                          'Switching Cost Reduction', 'Search Cost Reduction', 'No Info Campaign', 'Inc. All', 'Inc. All WH+']
        shares['Firm'] = ['ECS', 'EDF', 'Essent',
                          'ENINuon', 'Eneco', 'Lampiris', 'Other']
        # Setup table templates.
        table_row_header = ' & {val1} & {val2} & {val3} & {val4} & {val5} & {val6} & {val7} \\tabularnewline \n'
        table_row = '{firm} & {val1:.2f} & {val2:.2f} & {val3:.2f} & {val4:.2f} & {val5: .2f} & {val6: .2f} & {val7: .2f} \\tabularnewline \n'

        # Write to LaTeX table: Large table with observed vs. predicted
        with open(project_paths_join('OUT_TABLES', 'cf_shares_' + file_suffix+'.tex'), 'w') as tex_file:
            # Set a different types of table row with placeholders.
            # table_row_var = """\multirow{{2}}{{*}}{{{var}}} & {val1:.4f}{str1} & {val2:.2f} \\tabularnewline\n  & ({se1:.4f}) & \\tabularnewline\n"""
            # table_row_var_2 = """\multirow{{2}}{{*}}{{{var}}} & {val2:.4f}{str2} & {val3:.4f}{str3}  \\tabularnewline\n  & ({se2:.4f}) & ({se3:.4f}) \\tabularnewline\n"""
            # table_row_var_3 = """\multirow{{2}}{{*}}{{{var}}} &  &  & {val3:.4f}{str3}  \\tabularnewline\n  &  &  & ({se3:.4f}) \\tabularnewline\n"""
            # table_row_int = '{stat} & {val1:.4g} & {val2:.4g} & {val3:.4g}  \\tabularnewline\n'
            # table_row_float = '{stat} & {val1:.3f} & {val2:.3f} & {val3:.3f}  \\tabularnewline\n'

            # Top of table.
            tex_file.write('\\begin{tabular}{r|cc|cc|c|cc}\n\\toprule\n')
            # Header row with names and dependent variables.
            tex_file.write(table_row_header.format(val1=shares.columns.values[0],
                                                   val2=shares.columns.values[1],
                                                   val3=shares.columns.values[2],  val4=shares.columns.values[3],
                                                   val5=shares.columns.values[4],
                                                   val6=shares.columns.values[5],
                                                   val7=shares.columns.values[6]))
            tex_file.write('\\midrule \n')
            # Write shares to table.
            for index, row in shares.iterrows():
                tex_file.write(table_row.format(firm=shares['Firm'][index],
                                                val1=shares['Observed'][index],
                                                val2=shares['Predicted'][index],
                                                val3=shares['Switching Cost Reduction'][index],
                                                val4=shares['Search Cost Reduction'][index],
                                                val5=shares['No Info Campaign'][index],
                                                val6=shares['Inc. All'][index],
                                                val7=shares['Inc. All WH+'][index]))
            # Write bottom of table and notes.
            tex_file.write('\\midrule  \n')
            tex_file.write('\\multicolumn{8}{p{21.5cm}}{\\footnotesize{\\textit{Notes: The table summarizes observed (Column 1) and counterfactual (Columns 3 to 7) firm-level market shares averaged over all time periods. Column 1 and 2 compare market shares as observed in the data and predicted by our model. Column 3, 4, 5, 6, 7 display market shares from the switching costs reduction, the PCW search cost reduction, the no information campaign, the hypothetical monopolist offers all contracts at observed prices, the hypothetical monopolist offers all contracts at the indifference markup counterfactuals, respectively.}}}\\tabularnewline\n')

            # Bottom of table.
            tex_file.write('\\bottomrule\n\\end{tabular}\n')

        # Write to LaTeX table: Small table without observed versus predicted shares.
        with open(project_paths_join('OUT_TABLES', 'cf_shares_small_'+file_suffix+'.tex'), 'w') as tex_file:
            # Set a different types of table row with placeholders.
            table_row_header_small = ' & {val1} & {val2} & {val3} \\tabularnewline \n'
            table_row_small = '{firm} & {val1:.2f} & {val2:.2f} & {val3:.2f} \\tabularnewline\n'
            # Top of table.
            tex_file.write('\\begin{tabular}{r|c|cc}\n\\toprule\n')
            # Header row with names and dependent variables.
            tex_file.write(table_row_header_small.format(val1=shares.columns.values[0],
                                                         val2=shares.columns.values[2],  val3=shares.columns.values[3]))
            tex_file.write('\\midrule \n')
            # Write shares and HHI to table.
            for index, row in shares.iterrows():
                tex_file.write(table_row_small.format(firm=shares['Firm'][index],
                                                      val1=shares['Observed'][index],
                                                      val2=shares['Switching Cost Reduction'][index],
                                                      val3=shares['Search Cost Reduction'][index]
                                                      ))
            # Write bottom of table and notes.
            tex_file.write('\\midrule  \n')
            tex_file.write('\\multicolumn{4}{p{14.0cm}}{\\footnotesize{\\textit{Notes: The table summarizes observed (Column 1) and counterfactual (Columns 2 and 3) firm-level market shares averaged over all time periods. Column 2 and 3 display market shares from our switching costs reduction and the PCW search cost reduction counterfactuals,  respectively.}}}\\tabularnewline\n')
            # Bottom of table.
            tex_file.write('\\bottomrule\n\\end{tabular}\n')

        ################################################################################
        # 3. Plot indifference markups for different model specifications.
        # Read price data from csv-file.
        price_id_inc_1 = pd.read_csv(project_paths_join(
            'OUT_ANALYSIS', 'cf_inc_only1_price_id_' + file_suffix+'.csv'))
        price_id_inc_2 = pd.read_csv(project_paths_join(
            'OUT_ANALYSIS', 'cf_inc_only_2_price_id_' + file_suffix+'.csv'))
        price_id_inc_2_nosc = pd.read_csv(project_paths_join(
            'OUT_ANALYSIS', 'cf_inc_2_nosc_price_id_' + file_suffix+'.csv'))
        price_id_inc_all = pd.read_csv(project_paths_join(
            'OUT_ANALYSIS', 'cf_inc_all_price_id_' + file_suffix+'.csv'))
        # General plot prep.
        # Construct data axis.
        date_num = mdates.datestr2num(['2012-02', '2012-03', '2012-04', '2012-05', '2012-06',
                                       '2012-07', '2012-08', '2012-09', '2012-10', '2012-11', '2012-12',
                                       '2013-01', '2013-02', '2013-03', '2013-04', '2013-05', '2013-06',
                                       '2013-07', '2013-08', '2013-09', '2013-10', '2013-11', '2013-12',
                                       '2014-01', '2014-02', '2014-03', '2014-04', '2014-05', '2014-06',
                                       '2014-07', '2014-08', '2014-09', '2014-10', '2014-11', '2014-12',
                                       '2015-01', '2015-02', '2015-03', '2015-04', '2015-05', '2015-06',
                                       '2015-07', '2015-08', '2015-09', '2015-10', '2015-11', '2015-12',
                                       '2016-01', '2016-02', '2016-03', '2016-04', '2016-05', '2016-06'])

        # Define the basic characteristics of the figure.
        pyplot.clf()
        fig, ax = pyplot.subplots(figsize=(8, 6))
        ax = pyplot.subplot(111)
        pyplot.grid(True)
        pyplot.title(
            'Evolution of observed price and hypothetical indifference price path')
        # Plot price index datafor each operator.
        pyplot.plot_date(date_num, price_id_inc_1['Indiff. Markup'], 'r--',
                         label='Inc. only conv.', xdate=True, ydate=False)
        pyplot.plot_date(date_num, price_id_inc_2['Indiff. Markup'], 'b--',
                         label='Inc. conv. & green', xdate=True, ydate=False)
        pyplot.plot_date(date_num, price_id_inc_2_nosc['Indiff. Markup'], 'g--',
                         label='Inc. conv. & green, no SC', xdate=True, ydate=False)
        pyplot.plot_date(date_num, price_id_inc_all['Indiff. Markup'], 'y--',
                         label='Inc. all', xdate=True, ydate=False)
        pyplot.plot_date(date_num, price_id_inc_1['Observed'], 'k-',
                         label='Observed', xdate=True, ydate=False)

        ax.set_ylabel('Monthly Price (in EUR)')
        # ax.set_ylim([0,1])
        pyplot.xticks(np.arange(min(date_num), max(date_num)+1, 350))
        # Shink current axis by 17%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + 0.13, box.width, box.height * 0.83])
        ax.set_xlim([mdates.datestr2num('2012-02'),
                     mdates.datestr2num('2016-06')])
        # pyplot.xticks(np.arange(min(date_num), max(date_num)+180, 365.25))
        ax.xaxis.set_major_locator(mdates.MonthLocator(
            interval=10))  # to get a tick every 15 minutes
        ax.xaxis.set_major_formatter(
            mdates.DateFormatter('%b-%Y'))  # optional formatting
        # Put a legend to the right of the current axis
        # pyplot.legend(loc = 'upper center', ncol=2)
        ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.28), ncol=3)
        # Save the figure.
        pyplot.savefig(project_paths_join(
            'OUT_FIGURES', 'cf_price_id_' + file_suffix + '.pdf'), bbox_inches='tight')


        # Print average markups for different cases.
        print(markup_print.format(mt=model_type, gmm=gmm_type,
                                  ioc=100*id_mu_data['id_mu_avg']['Inc only conv.'],
                                  io=100*id_mu_data['id_mu_avg']['Inc only'],
                                  ionosc=100*id_mu_data['id_mu_avg']['Inc only (no SC)'], 
                                  iolcnosc=100*id_mu_data['id_mu_avg']['Inc only (LC, no SC)'], 
                                  incallnosc=100*id_mu_data['id_mu_avg']['Inc all'],
                                  obs=100*id_mu_data['id_mu_avg']['Obs']))
