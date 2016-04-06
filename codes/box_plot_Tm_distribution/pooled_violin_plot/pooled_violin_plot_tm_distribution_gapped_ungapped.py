#!/usr/bin/python

### plot amplicons for various size ranges
### Updated: version-1.03 02/19/2016
### Property of Wisser Lab at University of Delaware
### Author: Felix Francis (felixfrancier@gmail.com)


import matplotlib
import matplotlib.pyplot as plt
#matplotlib.use('Agg')  ### if running on a server and need to save the fig
import numpy as np
import pandas as pd
import statsmodels.api as sm
from matplotlib.ticker import NullFormatter



# input_file  = 'Tm_values.txt'
# input_file  = 'Tm_values_nogap.txt'
input_file_ungapped = 'Tm_values_UNgapped_only.txt'
input_file_gapped   = 'Tm_values_gapped_only.txt'


fields = ['Local_alignment_Tm', 'Gap_adjusted_end_filling_Tm']

df_ungapped = pd.read_table(input_file_ungapped , skipinitialspace=True, usecols=fields)      # skipinitialspace removes the spaces in the header. So ' Match_only_NN_Tm' becomes 'Match_only_NN_Tm'
df_gapped   = pd.read_table(input_file_gapped , skipinitialspace=True, usecols=fields)      # skipinitialspace removes the spaces in the header. So ' Match_only_NN_Tm' becomes 'Match_only_NN_Tm'



ungapped_local      = df_ungapped['Local_alignment_Tm'].values.tolist()
ungapped_end_filled  = df_ungapped['Gap_adjusted_end_filling_Tm'].values.tolist()

gapped_local      = df_gapped['Local_alignment_Tm'].values.tolist()
gapped_end_filled  = df_gapped['Gap_adjusted_end_filling_Tm'].values.tolist()





'''
#################################
### Tm DIFFERENCE VIOLIN PLOT
#################################

ungapped_difference = [a - b for a, b in zip(ungapped_local, ungapped_end_filled)]
gapped_difference   = [a - b for a, b in zip(gapped_local, gapped_end_filled)]

list = [ungapped_difference, gapped_difference ]

#fig = plt.figure()
fig = plt.figure(figsize=(8, 8), dpi=500)
ax = fig.add_subplot(111)
#sm.graphics.violinplot(list, ax=ax, plot_opts={'violin_fc':'grey', 'violin_alpha':0.25, 'violin_width':1, 'cutoff_val': 2})
sm.graphics.violinplot(list, ax=ax, positions=[.2,.55], plot_opts={'violin_fc':'grey', 'violin_alpha':0.25, 'violin_width':1})


min_list = []
max_list = []
for i in list:
    min_val = min(i)
    max_val = max(i)
    min_list.append(min_val)
    max_list.append(max_val)    
min_all     = min(min_list) -10
max_all     = max(max_list) +10



plt.xticks([.2,.55], ['Ungapped', 'Gapped'])
#plt.title("Comparison of Tm from Local alignment vs End filling")
plt.ylim([min_all, max_all])
plt.plot((.375, .375), (min_all, max_all), 'k-')
plt.xlim(0,.75)
plt.xlabel("HSE true hybridization prediction effect on Tm")
plt.ylabel("Tm difference")

#plt.show()
fig.savefig('difference_violin.png', dpi=500)
'''





#################################
### ALL 4 VIOLIN PLOT
#################################

list = [ungapped_local, ungapped_end_filled, gapped_local, gapped_end_filled ]
#
#print min(ungapped_local), min(ungapped_end_filled), min(gapped_local), min(gapped_end_filled)

min_list = []
max_list = []
for i in list:
    min_val = min(i)
    max_val = max(i)
    min_list.append(min_val)
    max_list.append(max_val)    
min_all     = min(min_list) -10
max_all     = max(max_list) +10



#fig = plt.figure()
fig = plt.figure(figsize=(13, 9), dpi=500)
ax = fig.add_subplot(111)
sm.graphics.violinplot(list, ax=ax, plot_opts={'violin_fc':'grey', 'violin_alpha':0.25, 'violin_width':1, 'cutoff_val': 2})
#sm.graphics.violinplot(list, ax=ax, plot_opts={'violin_fc':'grey', 'violin_alpha':0.25, 'violin_width':1})

plt.xticks([1, 2, 3, 4], ['Ungapped local', 'Ungapped end-filled', 'Gapped local', 'Gap adjusted & end-filled'])
#plt.title("Comparison of Tm from Local alignment vs End filling")
plt.ylim([min_all, max_all])
#plt.xlim(.8,2.2)
#plt.xlabel("HSE Alignment adjust method comparsion")
plt.ylabel("Tm")

plt.plot((2.5, 2.5), (min_all, max_all), 'k-')

#plt.show()

fig.savefig('all_four_violin.png', dpi=500)



