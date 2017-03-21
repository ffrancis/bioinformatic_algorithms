#!/usr/bin/env python

'''

Description:    Pandas vectorized solution to efficiently correct genotype calls based on coverage.
		If coverage < a user defined 'threshold', then the corresponding genotype values are changed to 0.


'''

#Time to run the code: start timer
import time
t0 = time.time()

import pandas
import numpy as np

############################################################
### INPUTS
############################################################

#input_file	= '/reads_to_genotype_data_small2.csv'
input_file	= '../data/reads_to_genotype_data_small.csv'
#input_file	= '../data/reads_to_genotype_data.csv'
threshold	= 12
output		= 'output.csv'



############################################################
### RUN
############################################################
df = pandas.read_csv(input_file)

n_samp	= (len(list(df.columns.values)) -5 )/2

col_names = list(df.columns.values)
col_names = col_names[5:]

df = df.replace(['.'], 0) 
#df = df.applymap(str)
df = pandas.to_numeric(df, errors='raise') 

#read_depth, genotype = 'HM160216_01_CML258.DP', 'HM160216_01_CML258.GT'
#df1.loc[df1[read_depth] < threshold, genotype] = 0
#
#print df1

for i in xrange(n_samp):
	read_depth, genotype = col_names[i], col_names[i+n_samp]
	df.loc[df[read_depth] < threshold, genotype] = 0

df.to_csv(output, sep=',')






#Time to run the code: end timer
t1 = time.time()
total = t1-t0
print (total)




