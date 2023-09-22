"""

Plot relation between supplier churn and vtest usage.

"""
import numpy as np
from bld.project_paths import project_paths_join
from matplotlib import pyplot
import matplotlib.dates as mdates
import datetime as dt
import pandas as pd
# Set up basic log file.
import logging
logging.basicConfig(filename=project_paths_join('OUT_ANALYSIS','awarenessswitching.log'),filemode='w',level=logging.INFO)
# Read price data from csv-file.
as_data = pd.read_csv(project_paths_join('OUT_DATA', 'awareness_switching_data.csv'))


# Construct data axis.
date_str = mdates.datestr2num(['2012-01','2012-02','2012-03','2012-04','2012-05','2012-06', \
                               '2012-07','2012-08','2012-09','2012-10','2012-11','2012-12', \
                               '2013-01','2013-02','2013-03','2013-04','2013-05','2013-06', \
                               '2013-07','2013-08','2013-09','2013-10','2013-11','2013-12', \
                               '2014-01','2014-02','2014-03','2014-04','2014-05','2014-06', \
                               '2014-07','2014-08','2014-09','2014-10','2014-11','2014-12', \
                               '2015-01','2015-02','2015-03','2015-04','2015-05','2015-06', \
                               '2015-07','2015-08','2015-09','2015-10','2015-11','2015-12', \
                               '2016-01','2016-02','2016-03','2016-04','2016-05','2016-06'])
date_str = ['2012-01','2012-02','2012-03','2012-04','2012-05','2012-06', \
                               '2012-07','2012-08','2012-09','2012-10','2012-11','2012-12', \
                               '2013-01','2013-02','2013-03','2013-04','2013-05','2013-06', \
                               '2013-07','2013-08','2013-09','2013-10','2013-11','2013-12', \
                               '2014-01','2014-02','2014-03','2014-04','2014-05','2014-06', \
                               '2014-07','2014-08','2014-09','2014-10','2014-11','2014-12', \
                               '2015-01','2015-02','2015-03','2015-04','2015-05','2015-06', \
                               '2015-07','2015-08','2015-09','2015-10','2015-11','2015-12', \
                               '2016-01','2016-02','2016-03','2016-04','2016-05','2016-06']

date_num = [dt.datetime.strptime(d,'%Y-%m').date() + dt.timedelta(days=30.475) for d in date_str]


# Define the basic characteristics of the figure.
pyplot.clf()
pyplot.plot()
pyplot.grid(True)
#pyplot.ylabel('Market shares in %')
# pyplot.title('Supplier Market Shares - Observed Data')

fig, ax1 = pyplot.subplots()
pyplot.grid(True)
ax2 = ax1.twinx()

# Plot information indicator on right axis.
ax1.plot(date_num, 100 *as_data['vtest'],'b-',label='Information indicator', linewidth=1.2)
ax1.set_ylim([0,13])
ax1.set_xlim([dt.date(2012, 1, 1),dt.date(2016, 6, 1)])

ax1.set_ylabel('Share of fully informed consumers (in %)')
# Plot switching rates on left axis.
ax2.plot(date_num, 100 * as_data['switching'], 'r-',label='Supplier churn rate', linewidth=0.8)
ax2.set_ylim([0,12])
ax2.set_ylabel('Total monthly churn rate (in %)')
ax2.set_xlim([dt.date(2012, 1, 1),dt.date(2016, 6, 1)])


myFmt = mdates.DateFormatter('%Y')
ax1.xaxis.set_major_formatter(myFmt)
start, end = ax1.get_xlim()
# ax1.xaxis.set_ticks(np.arange(start, end, 365))
n = 2  # Keeps every 7th label
[l.set_visible(False) for (i,l) in enumerate(ax1.xaxis.get_ticklabels()) if i % n != 0]
# pyplot.xticks(np.arange(min(date_num), max(date_num)+1, 350))

# Shink current axis by 17%
box = ax1.get_position()
ax1.set_position([box.x0, box.y0 + 0.13 , box.width * 0.95, box.height* 0.83])

# Put a legend to the right of the current axis
# ax.legend(loc='lower center', bbox_to_anchor=(1, 0.5))
ax1.legend(loc='lower center', bbox_to_anchor=(0.2, -0.28),ncol=1)
ax2.legend(loc='lower center', bbox_to_anchor=(0.8, -0.28),ncol=1)

# # Save the figure.
pyplot.savefig(project_paths_join('OUT_FIGURES', 'plotawarenessswitching.pdf'), bbox_inches ='tight')

## Look at raw correlation between PCW usage and aggregate supplier churn.
# Create year indicator in data.
as_data['year_str'] = as_data['month'].astype(str).str[0:4]
# Compute lagged values of churn rate and vtest usage.
# Probably lagged churn rate does nto make sense.
as_data['vtest_lag1'] = as_data['vtest'].shift(1)
as_data['vtest_lag2'] = as_data['vtest'].shift(2)
as_data['vtest_lag3'] = as_data['vtest'].shift(3)
# Compute moving averages of variables.
as_data['churn_ma'] = as_data['switching'].rolling(window=4).mean()
as_data['vtest_ma'] = as_data['vtest'].rolling(window=4).mean()
# Lage moving average vtest usage.
as_data['vtest_lag1_ma'] = as_data['vtest_ma'].shift(1)
as_data['vtest_lag2_ma'] = as_data['vtest_ma'].shift(2)
as_data['vtest_lag3_ma'] = as_data['vtest_ma'].shift(3)

# Compute correlations.
# Lagged PCW usage.
corr_0v = as_data['switching'].corr(as_data['vtest'])
corr_1v = as_data['switching'].corr(as_data['vtest_lag1'])
corr_2v = as_data['switching'].corr(as_data['vtest_lag2'])
corr_3v = as_data['switching'].corr(as_data['vtest_lag3'])

# Use moving average variables.
corr_0v = as_data['churn_ma'].corr(as_data['vtest_ma'])
corr_1v = as_data['churn_ma'].corr(as_data['vtest_lag1_ma'])
corr_2v = as_data['churn_ma'].corr(as_data['vtest_lag2_ma'])
corr_3v = as_data['churn_ma'].corr(as_data['vtest_lag3_ma'])


# logging.info(':Correlation when lagged switching is used:')
# logging.info('Contemporaneous correlation:')
# logging.info(corr_0s)
# logging.info('Between contemporaneous vtest and 1-month lagged switching:')
# logging.info(corr_1s)
# logging.info('Between contemporaneous vtest and 2-month lagged switching:')
# logging.info(corr_2s)
# logging.info('\n\n')
logging.info(':CORRELATION AVERAGED ACROSS ALL YEARS:')
logging.info('Contemporaneous correlation:')
logging.info(corr_0v)
logging.info('Between contemporaneous switching and 1-month lagged vtest:')
logging.info(corr_1v)
logging.info('Between contemporaneous switching and 2-month lagged vtest:')
logging.info(corr_2v)
logging.info('Between contemporaneous switching and 3-month lagged vtest:')
logging.info(corr_3v)


# Look at correlation separately for different years.
year_idx = ['2012','2013','2014', '2015', '2016']
for i, var in enumerate(year_idx):
    logging.info('\n\n')
    logging.info('Computing correlation for year:')
    logging.info(var)
    logging.info('\n')
    data_y = as_data[as_data.year_str==var]
    # Lagged switching.
#     corr_0sy = data_y['switching'].corr(data_y['vtest'])
#     corr_1sy = data_y['switching_lag1'].corr(data_y['vtest'])
#     corr_2sy = data_y['switching_lag2'].corr(data_y['vtest'])
    # Lagged PCW usage.
    corr_0vy = data_y['churn_ma'].corr(data_y['vtest_ma'])
    corr_1vy = data_y['churn_ma'].corr(data_y['vtest_lag1_ma'])
    corr_2vy = data_y['churn_ma'].corr(data_y['vtest_lag2_ma'])
    corr_3vy = data_y['churn_ma'].corr(data_y['vtest_lag3_ma'])

#     logging.info(':Correlation when lagged switching is used:')
#     logging.info('Contemporaneous correlation:')
#     logging.info(corr_0sy)
#     logging.info('Between contemporaneous vtest and 1-month lagged switching:')
#     logging.info(corr_1sy)
#     logging.info('Between contemporaneous vtest and 2-month lagged switching:')
#     logging.info(corr_2sy)
#     logging.info('\n')
#     logging.info(':Correlation when lagged vtest is used:')
    logging.info('Contemporaneous correlation:')
    logging.info(corr_0vy)
    logging.info('Between contemporaneous switching and 1-month lagged vtest:')
    logging.info(corr_1vy)
    logging.info('Between contemporaneous switching and 2-month lagged vtest:')
    logging.info(corr_2vy)
    logging.info('Between contemporaneous switching and 3-month lagged vtest:')
    logging.info(corr_3vy)