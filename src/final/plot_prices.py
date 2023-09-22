"""

Plot evolution of conventional and green prices.

"""

import numpy as np
from bld.project_paths import project_paths_join
from matplotlib import pyplot
import matplotlib.dates as mdates
import datetime as dt
import pandas as pd

# Read price data from csv-file.
cp_data = pd.read_csv(project_paths_join('OUT_DATA', 'prices_data.csv'))
cp_data_conv = cp_data[cp_data['contract_type'] == '1. conventional']
cp_data_green = cp_data[cp_data['contract_type'] == '3. green']

# Extract data for each supplier (conventional contracts).
cp_ecs_conv = cp_data_conv[cp_data_conv['firm'] == 'ECS']
cp_edf_conv = cp_data_conv[cp_data_conv['firm'] == 'EDF']
cp_eneco_conv = cp_data_conv[cp_data_conv['firm'] == 'Eneco']
cp_eni_conv = cp_data_conv[cp_data_conv['firm'] == 'Eni']
cp_essent_conv = cp_data_conv[cp_data_conv['firm'] == 'Essent']
cp_lampiris_conv = cp_data_conv[cp_data_conv['firm'] == 'Lampiris']
cp_other_conv = cp_data_conv[cp_data_conv['firm'] == 'Other']

# Extract data for each supplier (green contracts).
cp_ecs_green = cp_data_green[cp_data_green['firm'] == 'ECS']
cp_edf_green = cp_data_green[cp_data_green['firm'] == 'EDF']
cp_eneco_green = cp_data_green[cp_data_green['firm'] == 'Eneco']
cp_eni_green = cp_data_green[cp_data_green['firm'] == 'Eni']
cp_essent_green = cp_data_green[cp_data_green['firm'] == 'Essent']
cp_lampiris_green = cp_data_green[cp_data_green['firm'] == 'Lampiris']
cp_other_green = cp_data_green[cp_data_green['firm'] == 'Other']

# Construct data axis.
date_str = ['2012-01', '2012-02', '2012-03', '2012-04', '2012-05', '2012-06',
            '2012-07', '2012-08', '2012-09', '2012-10', '2012-11', '2012-12',
            '2013-01', '2013-02', '2013-03', '2013-04', '2013-05', '2013-06',
            '2013-07', '2013-08', '2013-09', '2013-10', '2013-11', '2013-12',
            '2014-01', '2014-02', '2014-03', '2014-04', '2014-05', '2014-06',
            '2014-07', '2014-08', '2014-09', '2014-10', '2014-11', '2014-12',
            '2015-01', '2015-02', '2015-03', '2015-04', '2015-05', '2015-06',
            '2015-07', '2015-08', '2015-09', '2015-10', '2015-11', '2015-12',
            '2016-01', '2016-02', '2016-03', '2016-04', '2016-05', '2016-06']
date_num = [dt.datetime.strptime(
    d, '%Y-%m').date() + dt.timedelta(days=30.475) for d in date_str]

# Plot conventional prices.
# Define the basic characteristics of the figure.
pyplot.clf()
pyplot.plot()
fig, ax = pyplot.subplots()
pyplot.grid(True)
# Plot information indicator on right axis.
ax.set_title('Conventional contracts')
ax.plot(date_num, cp_ecs_conv['price'], 'b-', label='ECS', linewidth=1.2)
ax.plot(date_num, cp_edf_conv['price'], 'r-', label='EDF', linewidth=1.0)
# ax.plot(date_num, cp_eneco_conv['price'],'r-',label='Eneco', linewidth=1.0)
ax.plot(date_num, cp_eni_conv['price'], 'g-', label='Eni', linewidth=1.0)
ax.plot(date_num, cp_essent_conv['price'], 'y-', label='Essent', linewidth=1.0)
# ax.plot(date_num, cp_lampiris_conv['price'],'r-',label='Lampiris', linewidth=1.0)
ax.plot(date_num, cp_other_conv['price'], 'k--', label='Other', linewidth=0.7)
ax.set_ylim([20, 40])
ax.set_xlim([dt.date(2012, 1, 1), dt.date(2016, 6, 1)])
ax.set_ylabel('Average real contract price (in EUR)')
myFmt = mdates.DateFormatter('%Y')
ax.xaxis.set_major_formatter(myFmt)
start, end = ax.get_xlim()
# Keeps every 2nd label
n = 2
[l.set_visible(False) for (i, l) in enumerate(
    ax.xaxis.get_ticklabels()) if i % n != 0]
# Shink current axis by 17%
box = ax.get_position()
ax.set_position([box.x0, box.y0 + 0.13, box.width, box.height * 0.83])
# Put a legend to the right of the current axis
ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.35), ncol=3)
# Save the figure.
pyplot.savefig(project_paths_join(
    'OUT_FIGURES', 'plotconvprices.pdf'), bbox_inches='tight')

# Plot green prices.
# Define the basic characteristics of the figure.
pyplot.clf()
pyplot.plot()
fig, ax = pyplot.subplots()
pyplot.grid(True)
# Plot information indicator on right axis.
ax.set_title('Green contracts')
ax.plot(date_num, cp_ecs_green['price'], 'b-', label='ECS', linewidth=1.2)
ax.plot(date_num, cp_edf_green['price'], 'r-', label='EDF', linewidth=1.0)
ax.plot(date_num, cp_eneco_green['price'], 'm-', label='Eneco', linewidth=1.0)
ax.plot(date_num, cp_eni_green['price'], 'g-', label='Eni', linewidth=1.0)
ax.plot(date_num, cp_essent_green['price'],
        'y-', label='Essent', linewidth=1.0)
ax.plot(date_num, cp_lampiris_green['price'],
        'c-', label='Lampiris', linewidth=1.0)
# ax.plot(date_num, cp_other_green['price'],'k--',label='Other', linewidth=0.7)
ax.set_ylim([20, 40])
ax.set_xlim([dt.date(2012, 1, 1), dt.date(2016, 6, 1)])
ax.set_ylabel('Average real contract price (in EUR)')
myFmt = mdates.DateFormatter('%Y')
ax.xaxis.set_major_formatter(myFmt)
start, end = ax.get_xlim()
# Keeps every 2nd label
n = 2
[l.set_visible(False) for (i, l) in enumerate(
    ax.xaxis.get_ticklabels()) if i % n != 0]
# Shink current axis by 17%
box = ax.get_position()
ax.set_position([box.x0, box.y0 + 0.13, box.width, box.height * 0.83])
# Put a legend to the right of the current axis
ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.35), ncol=3)
# Save the figure.
pyplot.savefig(project_paths_join(
    'OUT_FIGURES', 'plotgreenprices.pdf'), bbox_inches='tight')


# Plot comparison of average retail prices and wholesale prices to illustrate distribution of markups in our application.
data_full = pd.read_csv(project_paths_join('OUT_DATA', 'new_master_data.csv'))
wh_price_corr = pd.read_csv(project_paths_join('OUT_DATA', 'wholesale_corr_export.csv'))
# Aggregate data to month level (average across contracts and firms)
data_agg = data_full.groupby('month').aggregate(lambda x: np.mean(x))
# Extract data only for incumbent conventional contract.
data_inc = data_full.loc[(data_full['firm'] == 1) &
                         (data_full['contract_type'] == 1)]
data_inc.set_index(data_inc['month'], drop=True, inplace=True)
# Compute average retail price (weighted by contract market share).
data_avg_rp = data_full[['month', 'mshare_contract', 'price']]
# Compute weighted mean of retail price (across firms and contracts).


def wm(x): return np.average(
    x, weights=data_avg_rp.loc[x.index, "mshare_contract"])


avg_rp = data_avg_rp.groupby('month').agg(mshare_contract=(
    "mshare_contract", "sum"), price_weighted_mean=("price", wm))
avg_rp_uw = data_avg_rp.groupby('month').agg(mshare_contract=(
    "mshare_contract", "sum"), price_unweighted_mean=("price", "mean"))
# Write to aggregate data frame.
#data_agg['price_wm'] = avg_rp['price_weighted_mean']
data_agg['price_wm'] = avg_rp_uw['price_unweighted_mean']
data_agg['price_inc'] = data_inc['price']

wh_price_corr['month'] = data_agg.index
wh_price_corr.set_index(wh_price_corr['month'], drop=True, inplace=True)

# Transform wholesale prices into correct unit as described in MATLAB code.
# Retail price is monthly expenditure for a yearly 3500kwh consumption.
# Wholesale prices are for one mwh of electricity.
# Create vector for VAT.
T = 54
vat_rate = 1.21 * np.ones(T)
vat_rate[27:44] = 1.06
data_agg['wp_mwh'] = 1.0 * data_agg['wholesale_contract'].copy()
data_agg['wholesale_spot'] = vat_rate * \
    (data_agg['wholesale_spot'] * 3500 / 1000 / 12)
data_agg['wholesale_contract'] = vat_rate * \
    (data_agg['wholesale_contract'] * 3500 / 1000 / 12)
# Keep only columns that are needed in aggregate data.
data_agg = data_agg[['price_wm', 'price_inc',
                     'wholesale_spot', 'wholesale_contract', 'wp_mwh']]
data_agg['wh_peak'] = wh_price_corr['Peak']
data_agg['wh_late'] = wh_price_corr['Late']
data_agg['wh_early'] = wh_price_corr['Early']
# Compute markups: retail over wholesale price.
# Markups relative to wholesale price.
# data_agg['markup_inc'] = 100 * \
#     (data_agg['price_inc'] / data_agg['wholesale_spot'] - 1)
# data_agg['markup_spot'] = 100 * \
#     (data_agg['price_wm'] / data_agg['wholesale_spot'] - 1)
# data_agg['markup_contract'] = 100 * \
#     (data_agg['price_wm'] / data_agg['wholesale_contract'] - 1)
# Markups relative to retail price.
data_agg['markup_inc'] = 100 * \
    (data_agg['price_inc'] - data_agg['wh_peak']) / \
    data_agg['price_inc']
data_agg['markup_spot'] = 100 * \
    (data_agg['price_wm'] - data_agg['wh_peak']) / \
    data_agg['price_wm']
data_agg['markup_contract'] = 100 * \
    (data_agg['price_wm'] - data_agg['wholesale_contract']) / \
    data_agg['price_wm']


# Sanity check compute average retail price for one mwh.
data_agg['rp_mwh'] = data_agg['price_wm'] / 3500 * 1000 * 12 / vat_rate


################################################################################
# For comparison with ACER report, report everything in prices per mwh
data_agg['muw_mwh'] = 100 * ((data_agg['rp_mwh'] / data_agg['wh_peak']) - 1)
data_agg['mur_mwh'] = 100 * (data_agg['rp_mwh'] -
                             data_agg['wh_peak']) / data_agg['rp_mwh']
data_agg['muw_mwh_ma'] = data_agg['muw_mwh'].rolling(window=3).mean()
data_agg['mur_mwh_ma'] = data_agg['mur_mwh'].rolling(window=3).mean()
data_agg['margin_mwh'] = data_agg['rp_mwh'] - data_agg['wp_mwh']
# Define the basic characteristics of the figure.
pyplot.clf()
pyplot.plot()
fig, ax = pyplot.subplots()
# twin object for two different y-axis on the sample plot
ax2 = ax.twinx()
pyplot.grid(True)
# Plot information indicator on right axis.
#ax.set_title('Evolution of prices and markups')
ax.plot(date_num, data_agg['rp_mwh'], 'b-',
        label='Average retail price', linewidth=1.2)
#ax.plot(date_num, data_agg['price_inc'], 'r-', label='Incumbent retail price', linewidth=1.0)
ax.plot(date_num, data_agg['wh_peak'], 'g-',
        label='Wholesale Spot', linewidth=1.0)
# ax.plot(date_num, data_agg['margin_mwh'], 'r.',
#         label='Margin', markersize=1.2)
#ax.plot(date_num, data_agg['wholesale_contract'], 'y-', label='Wholesale Future', linewidth=1.0)
#ax2.plot(date_num, data_agg['markup_inc'], 'r--', label='Markup Incumbent', linewidth=0.7)
ax2.plot(date_num, data_agg['muw_mwh_ma'], 'r--',
         label='Avg. Markup (3-month MA)', linewidth=0.7)
# ax2.plot(date_num, data_agg['mur_mwh_ma'], 'y--',
#          label='Avg. Markup RP (MA)', linewidth=0.7)
# #ax.plot(date_num, data_agg['markup_contract'], 'r--', label='Avg. Markup Future', linewidth=0.7)
ax.set_ylim([0, 120])
ax2.set_ylim([40, 150])
ax.set_xlim([dt.date(2012, 1, 1), dt.date(2016, 6, 1)])
ax.set_ylabel('Prices (in EUR/mwh)')
ax2.set_ylabel('Average Markup (in percent)')
myFmt = mdates.DateFormatter('%Y')
ax.xaxis.set_major_formatter(myFmt)
start, end = ax.get_xlim()
# Keeps every 2nd label
n = 2
[l.set_visible(False) for (i, l) in enumerate(
    ax.xaxis.get_ticklabels()) if i % n != 0]
# Shink current axis by 17%
box = ax.get_position()
ax.set_position([box.x0, box.y0 + 0.13, box.width, box.height * 0.83])
ax2.set_position([box.x0, box.y0 + 0.06, box.width, box.height * 0.83])
# Put a legend to the right of the current axis
ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.23), ncol=3)
ax2.legend(loc='lower center', bbox_to_anchor=(0.5, -0.35), ncol=2)
# Save the figure.
pyplot.savefig(project_paths_join(
    'OUT_FIGURES', 'plot_prices_markups_mwh.pdf'), bbox_inches='tight')
np.mean(data_agg)
