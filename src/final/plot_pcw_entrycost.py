"""
	Plot the evolution of consumers' PCW search/entry costs
      based on model estimates.

      Output figures (one separate file foreach combination of model_list and gmm_list):
      - pcw_search_costs_FILESUFFIX.pdf

      Most importantly, this file formats Figure 1 on the main text.

"""
import numpy as np
from bld.project_paths import project_paths_join
from matplotlib import pyplot
import matplotlib.dates as mdates
import pandas as pd
# Select model types & GMM type to format.
model_list = [1, 4, 5, 6]
gmm_list = [
    [1, 2],
    [1, 2],
    [1],
    [1]]
model_print = r'Plotting results for implied PCW search costs based on model {mt:d} and GMM stage {gmm:d}...'
# Label columns.
kappa_columns = ['Average', 'Non-senior', 'Senior']
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
# Loop over different model specifications.
for iter, model_type in enumerate(model_list):
    # print(iter)
    for gmm_type in gmm_list[iter]:
        print(model_print.format(mt=model_type, gmm=gmm_type))
        file_suffix = str(gmm_type) + '_mod' + str(model_type)
        file_suffix_modelonly = '1' + '_mod' + str(model_type)
        # Read PCW search cost data from csv-file.
        kappa = pd.read_csv(project_paths_join(
            'OUT_ANALYSIS', 'kappa_EUR_' + file_suffix + '.csv'), header=None)
        kappa.columns = kappa_columns
        # Compute means over time.
        kappa['mean_average'] = [
            np.mean(kappa['Average'])]*len(kappa['Average'])
        kappa['mean_young'] = [
            np.mean(kappa['Non-senior'])]*len(kappa['Non-senior'])
        kappa['mean_old'] = [np.mean(kappa['Senior'])]*len(kappa['Senior'])
        # Define the basic characteristics of the figure.
        pyplot.clf()
        fig, ax = pyplot.subplots(figsize=(8, 6))
        ax = pyplot.subplot(111)
        pyplot.grid(True)
        
        # Plot price index datafor each operator.
        if model_type == 4:
            pyplot.title('PCW Search Cost over Time - Heterogeneous Frictions')
            pyplot.plot_date(
                date_num, kappa['Non-senior'], 'b-', label='Non-seniors',  xdate=True, ydate=False)
            pyplot.plot_date(
                date_num, kappa['mean_young'], 'b--', label='Non-seniors',  xdate=True, ydate=False)
            pyplot.plot_date(
                date_num, kappa['Senior'], 'r-', label='Seniors', xdate=True, ydate=False)
            pyplot.plot_date(
                date_num, kappa['mean_old'], 'r--', label='Seniors', xdate=True, ydate=False)
            pyplot.plot_date(
                date_num, kappa['Average'], 'k-', label='Average', xdate=True, ydate=False)
            pyplot.plot_date(
                date_num, kappa['mean_average'], 'k--', label='Average', xdate=True, ydate=False)
        else:
            pyplot.title('PCW Search Cost over Time - Homogeneous Frictions')
            pyplot.plot_date(
                date_num, kappa['Average'], 'k-', label='Average', xdate=True, ydate=False)
            pyplot.plot_date(
                date_num, kappa['mean_average'], 'k--', label='Average (over time)', xdate=True, ydate=False)

            # ax.plot(x, y_mean, label='Mean', linestyle='--')
        ax.set_ylabel('PCW Search Cost (in EUR)')
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
        if model_type == 4:
            mylabels = ['Non-senior', 'Non-senior (average over time)', 'Senior',
                        'Senior (average over time)', 'Average (over consumers)', 'Average (over time)']
        else:
            mylabels = ['Homogeneous PCW search cost', 'Average (over time)']
        ax.legend(labels=mylabels, loc='lower center',
                  bbox_to_anchor=(0.5, -0.28), ncol=3)
        # Save the figure.
        pyplot.savefig(project_paths_join(
            'OUT_FIGURES', 'pcw_search_costs_' + file_suffix + '.pdf'), bbox_inches='tight')
