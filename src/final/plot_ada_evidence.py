"""
	Plot the evidence for plausibility of 
    vtest-internet-advertising assumptions.

"""
import numpy as np
from bld.project_paths import project_paths_join
print("BLD IMPORT SUCCESSFUL!")
from matplotlib import pyplot
import matplotlib.dates as mdates
import pandas as pd



col = pd.Series(data=np.random.random_sample((1500,))*100)
dfInit = {}
idx = pd.date_range('2007-01-01', periods=1500, freq='M')
for i in range(100):
     dfInit[i] = col
dfInit['idx'] = idx
df = pd.DataFrame(dfInit).set_index('idx')

# Read price data from csv-file.
ada_data = pd.read_csv(project_paths_join('OUT_ANALYSIS', 'ada_evidence_data.csv'))


## First graph: PCW usage vs. ad spending.
# Sort values according to x-axis variable.
ada_data.sort_values('ad_spending_mly',ascending=True,inplace=True)
# Define the basic characteristics of the figure.
pyplot.clf()
fig, ax = pyplot.subplots(figsize=(8,6))
ax = pyplot.subplot(111)
pyplot.grid(True)
#pyplot.title('PCW usage vs. ad spending')
# Plot price index data for each operator.
pyplot.plot(ada_data['ad_spending_mly'], ada_data['vtest'] , 'b.',label='Data')
pyplot.plot(ada_data['ad_spending_mly'], ada_data['vt_ads_fit'],'r-',label='LOWESS fit')
ax.set_ylabel('Monthly share of PCW users')
ax.set_xlabel('Aggregate monthly advertising expenditure (in 10 mio. EUR)')
# Shink current axis by 17%
box = ax.get_position()
ax.set_position([box.x0, box.y0 + 0.13 , box.width, box.height* 0.83])
# Put a legend to the right of the current axis
ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.28),ncol=3)
# Save the figure.
pyplot.savefig(project_paths_join('OUT_FIGURES', 'ada_vtest_ads.pdf'), bbox_inches ='tight')

## Second graph: PCW usage vs. internet.
# Sort values according to x-axis variable.
ada_data.sort_values('fixed_broadband',ascending=True,inplace=True)
# Define the basic characteristics of the figure.
pyplot.clf()
fig, ax = pyplot.subplots(figsize=(8,6))
ax = pyplot.subplot(111)
pyplot.grid(True)
#pyplot.title('PCW usage vs. internet penetration')
# Plot price index data for each operator.
pyplot.plot(ada_data['fixed_broadband'], ada_data['vtest'] , 'b.',label='Data')
pyplot.plot(ada_data['fixed_broadband'], ada_data['vt_int_fit'],'r-',label='LOWESS fit')
ax.set_ylabel('Monthly share of PCW users')
ax.set_xlabel('Broadband internet penetration rate')
# Shink current axis by 17%
box = ax.get_position()
ax.set_position([box.x0, box.y0 + 0.13 , box.width, box.height* 0.83])
# Put a legend to the right of the current axis
ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.28),ncol=3)
# Save the figure.
pyplot.savefig(project_paths_join('OUT_FIGURES', 'ada_vtest_int.pdf'), bbox_inches ='tight')