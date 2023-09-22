"""

Plot the evolution of firm-level market shares for both observed and counterfactual shares.


"""

import numpy as np
from bld.project_paths import project_paths_join
from matplotlib import pyplot
import matplotlib.dates as mdates
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

# Read price data from csv-file.
ms_obs = pd.read_csv(project_paths_join('OUT_ANALYSIS', 'ms_firm_obs.csv'),header=None)
# Extract data columns for each operator.
ms_obs_ecs = ms_obs[[0]]
ms_obs_edf = ms_obs[[1]]
ms_obs_eneco = ms_obs[[2]]
ms_obs_eninuon = ms_obs[[3]]
ms_obs_essent = ms_obs[[4]]
ms_obs_lampiris = ms_obs[[5]]
ms_obs_other = ms_obs[[6]]


# Construct data axis.
date_num = mdates.datestr2num(['2012-02','2012-03','2012-04','2012-05','2012-06', \
                               '2012-07','2012-08','2012-09','2012-10','2012-11','2012-12', \
                               '2013-01','2013-02','2013-03','2013-04','2013-05','2013-06', \
                               '2013-07','2013-08','2013-09','2013-10','2013-11','2013-12', \
                               '2014-01','2014-02','2014-03','2014-04','2014-05','2014-06', \
                               '2014-07','2014-08','2014-09','2014-10','2014-11','2014-12', \
                               '2015-01','2015-02','2015-03','2015-04','2015-05','2015-06', \
                               '2015-07','2015-08','2015-09','2015-10','2015-11','2015-12', \
                               '2016-01','2016-02','2016-03','2016-04','2016-05','2016-06'])

# Define the basic characteristics of the figure.
pyplot.clf()
pyplot.plot()
pyplot.grid(True)
pyplot.ylabel('Market shares')
#pyplot.title('Supplier Market Shares - Observed Data')
ax = pyplot.subplot(111)
# Plot price index datafor each operator.
pyplot.plot_date(date_num, np.array(ms_obs_ecs), 'r-',label='ECS', xdate=True,ydate=False)
pyplot.plot_date(date_num, np.array(ms_obs_edf),'y-',label='EDF',  xdate=True,ydate=False)
pyplot.plot_date(date_num, np.array(ms_obs_eneco), 'b-',label='Eneco', xdate=True,ydate=False)
pyplot.plot_date(date_num, np.array(ms_obs_essent), 'm-',label='Essent', xdate=True,ydate=False)
pyplot.plot_date(date_num, np.array(ms_obs_eninuon),'g-',label='EniNUON', xdate=True,ydate=False)
pyplot.plot_date(date_num, np.array(ms_obs_lampiris),'c-',label='Lampiris', xdate=True,ydate=False)
# pyplot.plot_date(date_num, ms_obs_other, 'k--',label='Other', xdate=True,ydate=False)

ax.set_ylim([0,0.55])
pyplot.xticks(np.arange(min(date_num), max(date_num)+1, 350))
ax.set_xlim([mdates.datestr2num('2012-02'),mdates.datestr2num('2016-06')])
# pyplot.xticks(np.arange(min(date_num), max(date_num)+180, 365.25))
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=10))   #to get a tick every 15 minutes
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b-%Y'))     #optional formatting 
# Shink current axis by 15%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.83, box.height])
# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# # Save the figure.
pyplot.savefig(project_paths_join('OUT_FIGURES', 'plotmsobs.pdf'), bbox_inches ='tight')

###########################################
## PLOT MARKET SHARES FOR COUNTERFACTUALS.#
###########################################
# Select model types & GMM type to format.
model_list = [1, 4, 5, 6]
gmm_list = [
    [2],
    [2],
    [1],
    [1]]
model_print = r'Plotting counterfactual market shares for model {mt:d} and GMM stage {gmm:d}...\n'
for iter, model_type in enumerate(model_list):
    #print(iter)
    for gmm_type in gmm_list[iter]:
        print(model_print.format(mt=model_type, gmm=gmm_type))
        file_suffix = str(gmm_type) + '_mod' + str(model_type)
        file_suffix_modelonly = '_mod' + str(model_type)
        # For switching subsidy.
        # Read price data from csv-file.
        ms_subsidy = pd.read_csv(project_paths_join('OUT_ANALYSIS', 'ms_firm_subsidy'+file_suffix_modelonly+'.csv'),header=None)
        # Extract data columns for each operator.
        ms_subsidy_ecs = ms_subsidy[[0]]
        ms_subsidy_edf = ms_subsidy[[1]]
        ms_subsidy_eneco = ms_subsidy[[2]]
        ms_subsidy_eninuon = ms_subsidy[[3]]
        ms_subsidy_essent = ms_subsidy[[4]]
        ms_subsidy_lampiris = ms_subsidy[[5]]
        ms_subsidy_other = ms_subsidy[[6]]

        # Define the basic characteristics of the figure.
        pyplot.clf()
        pyplot.plot()
        pyplot.grid(True)
        pyplot.ylabel('Market shares')
        #pyplot.title('Supplier Market Shares - Reduced Switching Costs ')
        ax = pyplot.subplot(111)
        # Plot price index datafor each operator.
        pyplot.plot_date(date_num, np.array(ms_subsidy_ecs), 'r-',label='ECS', xdate=True,ydate=False)
        pyplot.plot_date(date_num, np.array(ms_subsidy_edf),'y-',label='EDF',  xdate=True,ydate=False)
        pyplot.plot_date(date_num, np.array(ms_subsidy_eneco), 'b-',label='Eneco', xdate=True,ydate=False)
        pyplot.plot_date(date_num, np.array(ms_subsidy_essent), 'm-',label='Essent', xdate=True,ydate=False)
        pyplot.plot_date(date_num, np.array(ms_subsidy_eninuon),'g-',label='EniNUON', xdate=True,ydate=False)
        pyplot.plot_date(date_num, np.array(ms_subsidy_lampiris),'c-',label='Lampiris', xdate=True,ydate=False)
        # pyplot.plot_date(date_num, ms_subsidy_other, 'k--',label='Other', xdate=True,ydate=False)

        ax.set_ylim([0,0.55])
        pyplot.xticks(np.arange(min(date_num), max(date_num)+1, 350))
        ax.set_xlim([mdates.datestr2num('2012-02'),mdates.datestr2num('2016-06')])
        # pyplot.xticks(np.arange(min(date_num), max(date_num)+180, 365.25))
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=10))   #to get a tick every 15 minutes
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%b-%Y'))     #optional formatting 
        # Shink current axis by 15%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.83, box.height])
        # Put a legend to the right of the current axis
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        # Save the figure.
        pyplot.savefig(project_paths_join('OUT_FIGURES', 'plotmssubsidy_'+file_suffix+'.pdf'), bbox_inches ='tight')

        # For perfect information.
        # Read price data from csv-file.
        ms_perfect_info = pd.read_csv(project_paths_join('OUT_ANALYSIS', 'ms_firm_perfect_info' + file_suffix_modelonly + '.csv'),header=None)
        # Extract data columns for each operator.
        ms_perfect_info_ecs = ms_perfect_info[[0]]
        ms_perfect_info_edf = ms_perfect_info[[1]]
        ms_perfect_info_eneco = ms_perfect_info[[2]]
        ms_perfect_info_eninuon = ms_perfect_info[[3]]
        ms_perfect_info_essent = ms_perfect_info[[4]]
        ms_perfect_info_lampiris = ms_perfect_info[[5]]
        ms_perfect_info_other = ms_perfect_info[[6]]

        # Define the basic characteristics of the figure.
        pyplot.clf()
        pyplot.plot()
        pyplot.grid(True)
        pyplot.ylabel('Market shares')
        #pyplot.title('Supplier Market Shares - Reduced PCW Search Costs')
        ax = pyplot.subplot(111)
        # Plot price index datafor each operator.
        pyplot.plot_date(date_num, np.array(ms_perfect_info_ecs), 'r-',label='ECS', xdate=True,ydate=False)
        pyplot.plot_date(date_num, np.array(ms_perfect_info_edf),'y-',label='EDF',  xdate=True,ydate=False)
        pyplot.plot_date(date_num, np.array(ms_perfect_info_eneco), 'b-',label='Eneco', xdate=True,ydate=False)
        pyplot.plot_date(date_num, np.array(ms_perfect_info_essent), 'm-',label='Essent', xdate=True,ydate=False)
        pyplot.plot_date(date_num, np.array(ms_perfect_info_eninuon),'g-',label='EniNUON', xdate=True,ydate=False)
        pyplot.plot_date(date_num, np.array(ms_perfect_info_lampiris),'c-',label='Lampiris', xdate=True,ydate=False)
        # pyplot.plot_date(date_num, ms_perfect_info_other, 'k--',label='Other', xdate=True,ydate=False)

        ax.set_ylim([0,0.55])
        pyplot.xticks(np.arange(min(date_num), max(date_num)+1, 350))
        ax.set_xlim([mdates.datestr2num('2012-02'),mdates.datestr2num('2016-06')])
        # pyplot.xticks(np.arange(min(date_num), max(date_num)+180, 365.25))
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=10))   #to get a tick every 15 minutes
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%b-%Y'))     #optional formatting 
        # Shink current axis by 15%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.83, box.height])
        # Put a legend to the right of the current axis
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        # Save the figure.
        pyplot.savefig(project_paths_join('OUT_FIGURES', 'plotmsperfectinfo_'+file_suffix+'.pdf'), bbox_inches ='tight')

