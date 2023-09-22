"""
	Plot the evolution of consumers' PCW usage as observed/predicted by our data/model.
	Plot the evoluation of firm-level churn rates as observed/predicted by the model.

    OUTPUT FIGURES (separate file for each combination of model_list and gmm_list):
    - gof_pcw_churn_FILESUFFIX.pdf
    - gof_pcw_churn_combined_FILESUFFIX.pdf

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
model_print = r'Formatting estimation results for model {mt:d} and GMM stage {gmm:d}...\n'

for iter, model_type in enumerate(model_list):
    for gmm_type in gmm_list[iter]:
        print(model_print.format(mt=model_type, gmm=gmm_type))
        file_suffix = str(gmm_type) + '_mod' + str(model_type)
        file_suffix_modelonly = '1' + '_mod' + str(model_type)
        # Start actual code here.
        data_filename = 'gf_raw_data_' + \
            str(gmm_type) + '_mod' + str(model_type) + '.csv'
        # Read price data from csv-file.
        gof_data = pd.read_csv(project_paths_join(
            'OUT_ANALYSIS', data_filename))

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

        # Observed vs. predicted PCW usage and churn in two subplots.
        pyplot.clf()
        fig, ax = pyplot.subplots(figsize=(8, 6))
        ax = pyplot.subplot(211)
        pyplot.grid(True)
        pyplot.title('Aggregate PCW usage over time')
        # Plot price index datafor each operator.
        pyplot.plot_date(date_num[1:], gof_data['pcwObs'][1:],
                        'r--', label='Data', xdate=True, ydate=False)
        pyplot.plot_date(date_num[1:], gof_data['pcwPred'][:-1],
                        'b--', label='Predicted',  xdate=True, ydate=False)
        # pyplot.plot_date(date_num, pcw_comp['Switching Subsidy'], 'g-',label='Reduced switching cost', xdate=True,ydate=False)
        # pyplot.plot_date(date_num, pcw_comp['Search Cost Reduction'],'y-',label='Reduced PCW search cost', xdate=True,ydate=False)
        ax.set_ylabel('Monthly share of PCW users')
        ax.set_ylim([0, 0.15])
        pyplot.xticks(np.arange(min(date_num), max(date_num)+1, 350))
        # Shink current axis by 17%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + 0.15, box.width, box.height * 0.83])
        ax.set_xlim([mdates.datestr2num('2012-02'), mdates.datestr2num('2016-06')])
        # pyplot.xticks(np.arange(min(date_num), max(date_num)+180, 365.25))
        ax.xaxis.set_major_locator(mdates.MonthLocator(
            interval=10))  # to get a tick every 15 minutes
        ax.xaxis.set_major_formatter(
            mdates.DateFormatter('%b-%Y'))  # optional formatting
        # Put a legend to the right of the current axis
        pyplot.legend(loc='upper center', ncol=2)
        ax.legend(loc='upper right', bbox_to_anchor=(0.65, -0.14), ncol=2)

        ax = pyplot.subplot(212)
        pyplot.grid(True)
        pyplot.title('Aggregate churn rates over time')
        # Plot price index datafor each operator.
        pyplot.plot_date(date_num[1:], gof_data['churnObs'][1:],
                        'r--', label='Data', xdate=True, ydate=False)
        pyplot.plot_date(date_num[1:], gof_data['churnPred'][:-1],
                        'b--', label='Predicted',  xdate=True, ydate=False)
        # pyplot.plot_date(date_num, pcw_comp['Switching Subsidy'], 'g-',label='Reduced switching cost', xdate=True,ydate=False)
        # pyplot.plot_date(date_num, pcw_comp['Search Cost Reduction'],'y-',label='Reduced PCW search cost', xdate=True,ydate=False)
        ax.set_ylabel('Monthly firm-level churn rate')
        ax.set_ylim([0, 0.06])
        pyplot.xticks(np.arange(min(date_num), max(date_num)+1, 350))
        # Shink current axis by 17%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + 0.13, box.width, box.height * 0.83])
        ax.set_xlim([mdates.datestr2num('2012-02'), mdates.datestr2num('2016-06')])
        # pyplot.xticks(np.arange(min(date_num), max(date_num)+180, 365.25))
        ax.xaxis.set_major_locator(mdates.MonthLocator(
            interval=10))  # to get a tick every 15 minutes
        ax.xaxis.set_major_formatter(
            mdates.DateFormatter('%b-%Y'))  # optional formatting
        # Put a legend to the right of the current axis
        #pyplot.legend(loc='upper center', ncol=2)
        #ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.28), ncol=2)
        # Save the figure.
        pyplot.savefig(project_paths_join(
            'OUT_FIGURES', 'gof_pcw_churn_' + file_suffix +'.pdf'), bbox_inches='tight')

        # Observed vs. predicted PCW usage and churn in one combined plot.
        # Probably this is silly...so just not do it.
        pyplot.clf()
        fig, ax = pyplot.subplots(figsize=(8, 6))
        ax = pyplot.subplot(111)
        pyplot.grid(True)
        pyplot.title('Aggregate PCW usage and churn rates over time')
        # Plot price index datafor each operator.
        pyplot.plot_date(date_num[1:], gof_data['pcwObs'][1:],
                        'r-', label='Observed PCW usage', xdate=True, ydate=False)
        pyplot.plot_date(date_num[1:], gof_data['pcwPred'][:-1],
                        'r--', label='Predicted PCW usage',  xdate=True, ydate=False)
        # Plot price index datafor each operator.
        pyplot.plot_date(date_num[1:], gof_data['churnObs'][1:],
                        'b-', label='Observed churn', xdate=True, ydate=False)
        pyplot.plot_date(date_num[1:], gof_data['churnPred'][:-1],
                        'b--', label='Predicted churn',  xdate=True, ydate=False)

        ax.set_ylabel('Monthly share of PCW users and churn rate')
        ax.set_ylim([0, 0.15])
        pyplot.xticks(np.arange(min(date_num), max(date_num)+1, 350))
        # Shink current axis by 17%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + 0.15, box.width, box.height * 0.83])
        ax.set_xlim([mdates.datestr2num('2012-02'), mdates.datestr2num('2016-06')])
        # pyplot.xticks(np.arange(min(date_num), max(date_num)+180, 365.25))
        ax.xaxis.set_major_locator(mdates.MonthLocator(
            interval=10))  # to get a tick every 15 minutes
        ax.xaxis.set_major_formatter(
            mdates.DateFormatter('%b-%Y'))  # optional formatting
        # Put a legend to the right of the current axis
        pyplot.legend(loc='upper center', ncol=2)
        ax.legend(loc='upper right', bbox_to_anchor=(0.65, -0.14), ncol=2)
        # Save the figure.
        pyplot.savefig(project_paths_join(
            'OUT_FIGURES', 'gof_pcw_churn_combined_' + file_suffix + '.pdf'), bbox_inches='tight')
