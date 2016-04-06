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
input_file  = 'Tm_values_UNgapped_only.txt'





# fields = ['Match_only_NN_Tm', 'MisMatch_only_NN_Tm']
# fields = ['Local alignment Tm', 'End filling based Tm']
# fields = ['Local alignment Tm', 'Gap adjusted end filling based Tm']
fields = ['Local_alignment_Tm', 'Gap_adjusted_end_filling_Tm']

df = pd.read_table(input_file, skipinitialspace=True, usecols=fields)      # skipinitialspace removes the spaces in the header. So ' Match_only_NN_Tm' becomes 'Match_only_NN_Tm'





Local_alignment_Tm  = df['Local_alignment_Tm'].values.tolist()
Gap_adjusted_end_filling_Tm  = df['Gap_adjusted_end_filling_Tm'].values.tolist()



#################################
### SCATTER HISTOGRAM
#################################



nullfmt = NullFormatter()         # no labels


x   = Local_alignment_Tm
y   = Gap_adjusted_end_filling_Tm 

min_x =min(Local_alignment_Tm)
min_y =min(Gap_adjusted_end_filling_Tm)

min_Tm = min(min_x, min_y)
max_Tm = max(min_x, min_y)
       
#print min_Tm, max_Tm
         
# definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
bottom_h = left_h = left + width + 0.02

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.2, height]


# start with a rectangular Figure
#plt.figure(1, figsize=(12, 12))
plt.figure(1, figsize=(8, 8), dpi=500)

axScatter = plt.axes(rect_scatter)

plt.xlabel('Local alignment Tm')
plt.ylabel('HSE end filling Tm')

axHistx = plt.axes(rect_histx)
axHisty = plt.axes(rect_histy)

# no labels
axHistx.xaxis.set_major_formatter(nullfmt)
axHisty.yaxis.set_major_formatter(nullfmt)

# the scatter plot:
axScatter.scatter(x, y, s=10, facecolors='none', edgecolors='black', alpha=0.1)




# now determine nice limits by hand:
binwidth = 5
#binwidth = 0.25
xymax = np.max([np.max(np.fabs(x)), np.max(np.fabs(y))])
lim = (int(xymax/binwidth) + 1) * binwidth


axScatter.set_xlim((-lim, lim))
axScatter.set_ylim((-lim, lim))


#axScatter.set_xlim((min_Tm-10, max_Tm+10))
#axScatter.set_ylim((min_Tm-10, max_Tm+10))


bins = np.arange(-lim, lim + binwidth, binwidth)

#print bins



# histogram
axHistx.hist(x, bins=bins, color = 'grey')
axHisty.hist(y, bins=bins, orientation='horizontal',color = 'grey')


axHistx.set_xlim(axScatter.get_xlim())
axHistx.set_ylim(0,20000)

axHisty.set_xlim(0,20000)
axHisty.set_ylim(axScatter.get_ylim())

plt.setp(plt.xticks()[1], rotation=270)



#plt.show()
plt.savefig('ungapped_scatter.png', dpi=500)





