### repeat quantifciation
### Only 25 kb primers will be selected
### Only ones with Tm difference <10C will be counted.
### The center point of the 25 bp primer will be used as the coordinate...so 13th bp position.
### TIMESTAMP_HSE_out1_1.csv has the alignments for the + strand primers.
### TIMESTAMP_HSE_out3_1.csv has the coordinate information for these primers.



############################################################
# Time to run the code: start timer
############################################################
import time
t0 = time.time()


############################################################
# Import
############################################################
import numpy as np
import pandas as pd
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv
from itertools import izip
from parameters import *
import seaborn as sns
import copy

### FUNCTIONS



# FUNCTION TO CALCUALTE PERCENTAGE ID BETWEEN TWO SEQUENCES
def pid(seq1, seq2):
    if len(seq1) == len(seq2):
        score = 0
        seq1 = (seq1).upper()
        seq2 = (seq2).upper()
        for i in xrange(len(seq1)):
            if seq1[i] != seq2[i]:
                score = score+1
        pid = (len(seq1) - score) / float(len(seq1))
        pid = pid*100
    else:
        ##pid = "different_lengths!"
        pid = 0
    return pid    

# FUNCTION TO CALCULATE GC CONTENT
def gc_content(sequence):
    sequence = sequence.upper()
    gc_count =  sequence.count('G')+sequence.count('C')
    sequence_len = len(sequence)
    gc_cont    = (float(gc_count)/sequence_len)*100
    return round(gc_cont, 2)  
     
# ADD ZEROS TO MISSING COORDINATES (USES start_pos, stop_pos FROM PARAMETER.PY)
def fill_in_missing_coords(x_list, y_list):
    for i in xrange(start_pos, stop_pos+1):
        if i in x_list:
            pass
        else:
            x_list.append(i)
            y_list.append(0)

# ADD nan TO MISSING COORDINATES (USES start_pos, stop_pos FROM PARAMETER.PY)
def fill_in_missing_coords_nan(x_list, y_list):
    for i in xrange(start_pos, stop_pos+1):
        if i in x_list:
            pass
        else:
            x_list.append(i)
            y_list.append(np.nan)
            
            
            
            
### INPUT FILES AND PATHS

file_path = './'

pid_threshold = 70

time_stamp = 'TA_2016-07-05T19_03_49_824843'        # TEST DATA 4KB LOCUS
# time_stamp = 'TA_2016-07-12T17_04_51_600848'

f_primer_alignments_no_em = time_stamp+"_HSE_out1_1.csv"
f_primer_alignments_em = time_stamp+"_HSE_out1_1_em.csv"
primer_coord_files_no_em  = time_stamp+"_HSE_out3_1.csv"
primer_coord_files_em  = time_stamp+"_HSE_out3_1_em.csv"





### CODE

df_f_primers_no_em = pd.read_csv(file_path + f_primer_alignments_no_em)
df_f_primers_em = pd.read_csv(file_path + f_primer_alignments_em)
primer_names_no_em = pd.read_csv(file_path + primer_coord_files_no_em)
primer_names_em = pd.read_csv(file_path + primer_coord_files_em)

# MERGE EM AND NON-EM DATAFRAMES

df_f_primers_frames = [df_f_primers_no_em, df_f_primers_em]
primer_names_frames = [primer_names_no_em, primer_names_em]

df_f_primers = pd.concat(df_f_primers_frames)
primer_names = pd.concat(primer_names_frames)




# CREATE DICTIONARY OF PRIMER SEQUENCE AND PRIMER COORDINATE
primer_coord_dict = dict(zip(primer_names.Primer, primer_names.Genome_start))




# CREATE DICTIONARY OF PRIMER COORDINATE AND PRIMER GC
primer_gc_dict = dict(zip(primer_names.Primer, primer_names.Primer_GC,))
# print primer_gc_dict





# IF PRIMER LENGTH NOT EQUAL TO 25bp, REPLACE THAT PRIMER SEQUENCE WITH NaN
df_f_primers['Primer'] = np.where( df_f_primers['Primer'].str.len() == 25, df_f_primers['Primer'], np.nan)
# REMOVE ROWS WITH NaN AS PRIMER SEQUENCE (PRIMER LENGTHS != 25)
df_f_primers = df_f_primers.dropna(subset=['Primer'])


df_f_primers_tm = df_f_primers.copy(deep=True)


# IF Tm DIFFERENCE <10, RECORD THE VALUE ELSE NaN UNDER Tm DIFFERENCE COLUMN
df_f_primers_tm['Tm_difference'] = np.where( df_f_primers_tm['Primer_NN_Tm'] - df_f_primers_tm['Gap_adjusted_end_filling_Tm'] < 10, (df_f_primers_tm['Primer_NN_Tm']) - (df_f_primers_tm['Gap_adjusted_end_filling_Tm']), np.nan)
# REMOVE ROWS WITH NaN AS Tm DIFFERENCE
df_f_primers_tm = df_f_primers_tm.dropna(subset=['Tm_difference'])



# CALCULATE THE PID BETWEEN PRIMER SEQ AND THERMOALIGN SEQUENCE
df_f_primers['p_id'] = df_f_primers.apply(lambda row: pid(row['Primer'], row['3&5_prime_match_extend']), axis=1)



# IF PID NOT EQUAL TO THRESHOLD, REPLACE THAT PRIMER SEQUENCE WITH NaN
df_f_primers['p_id'] = np.where( df_f_primers['p_id'] >= pid_threshold, df_f_primers['p_id'], np.nan)
# REMOVE ROWS WITH NaN AS PRIMER SEQUENCE (PID LESS THAN pid_threshold)
df_f_primers = df_f_primers.dropna(subset=['p_id'])


# print df_f_primers


# TM BASED PLOT INFO
# COUNT THE NUMBER OF OCCURENCE OF UNIQUE PRIMER SEQUENCES IN PRIMER COLUMN. THIS INDICATES THE NUMBER OF NON TARGET HITS EACH PRIMER HAS
counts_tm = df_f_primers_tm['Primer'].value_counts()
# REPLACE THE PRIMER SEQUENCES WITH THE CORRESPONDING GENOME CORRDINATES FROM primer_coord_dict
index_coord_list = [primer_coord_dict[x] for x in counts_tm.index.values.tolist()]

x_list  = index_coord_list
y_list_unsorted = counts_tm.values.tolist()
 
original_tm_list = copy.deepcopy(y_list_unsorted )




fill_in_missing_coords(x_list, y_list_unsorted)
    
y_list = [x for (y,x) in sorted(zip(x_list, y_list_unsorted))]
x_list.sort()







# PID BASED PLOT INFO
# COUNT THE NUMBER OF OCCURENCE OF UNIQUE PRIMER SEQUENCES IN PRIMER COLUMN. THIS INDICATES THE NUMBER OF NON TARGET HITS EACH PRIMER HAS
counts_pid = df_f_primers['Primer'].value_counts()
# REPLACE THE PRIMER SEQUENCES WITH THE CORRESPONDING GENOME CORRDINATES FROM primer_coord_dict
index_coord_list_pid = [primer_coord_dict[x] for x in counts_pid.index.values.tolist()]

x_list_pid  = index_coord_list_pid
y_list_unsorted_pid = counts_pid.values.tolist()
 
 
 
fill_in_missing_coords(x_list_pid, y_list_unsorted_pid)
    
y_list_pid = [x for (y,x) in sorted(zip(x_list_pid, y_list_unsorted_pid))]
x_list_pid.sort()




# GC CONTENT INFO
# COUNT THE NUMBER OF OCCURENCE OF UNIQUE PRIMER SEQUENCES IN PRIMER COLUMN. THIS INDICATES THE NUMBER OF NON TARGET HITS EACH PRIMER HAS
counts_pid = df_f_primers['Primer'].value_counts()
# REPLACE THE PRIMER SEQUENCES WITH THE CORRESPONDING GENOME CORRDINATES FROM primer_coord_dict
index_coord_list_pid = [primer_coord_dict[x] for x in counts_pid.index.values.tolist()]
primer_list_pid = counts_pid.index.values.tolist()
# COUNT THE NUMBER OF OCCURENCE OF UNIQUE PRIMER SEQUENCES IN PRIMER COLUMN. THIS INDICATES THE NUMBER OF NON TARGET HITS EACH PRIMER HAS
counts_tm = df_f_primers_tm['Primer'].value_counts()
# REPLACE THE PRIMER SEQUENCES WITH THE CORRESPONDING GENOME CORRDINATES FROM primer_coord_dict
index_coord_list = [primer_coord_dict[x] for x in counts_tm.index.values.tolist()]
primer_list_tm = counts_tm.index.values.tolist()


merged_primer_list = primer_list_pid + list(set(primer_list_tm) - set(primer_list_pid))


gc_list = [primer_gc_dict[x] for x in merged_primer_list]


x_list_gc  = index_coord_list_pid  # x axis data; gc_list = y axis data
y_list_unsorted_gc = gc_list

original_gc_list = copy.deepcopy(y_list_unsorted_gc)

 
fill_in_missing_coords_nan(x_list_gc, y_list_unsorted_gc)
    
y_list_gc = [x for (y,x) in sorted(zip(x_list_gc, y_list_unsorted_gc))]
x_list_gc.sort()





plt.figure(figsize=(20,16))

#### Set the font dictionaries (for plot title and axis titles)
title_font = { 'size':'20', 'color':'black', 'weight':'normal',
              'verticalalignment':'bottom'} # Bottom vertical alignment for more space
axis_font = { 'size':'20'}


plt.rc('xtick', labelsize=15) 
plt.rc('ytick', labelsize=20)

plt.plot(x_list,y_list, color="black", label='Tm based repeats')
plt.fill_between(x_list, y_list, color="grey")


plt.plot(x_list_pid,y_list_pid, color="red", alpha=0.5, label= '>= ' +str(pid_threshold) + ' p_id based repeats')
plt.fill_between(x_list_pid,y_list_pid, color="red", alpha=0.3)
# plt.plot(x_list,y_list, color="black", label='Tm based repeats')
# plt.fill_between(x_list, y_list, color="grey")
plt.xlabel('Coordinate', **axis_font)
plt.ylabel('Repeat Count', **axis_font)
plt.suptitle('Repeat count distribution for qSLB_3_33490673_24001', **title_font)
plt.xticks(rotation=70)
plt.ticklabel_format(useOffset=False, style='plain')

plt.legend(loc='upper left')
plt.savefig('Tm_pid_hits_plot.png', dpi=300)


############ TWO Y AXIS ##############

# plt = plt.figure(figsize=(20,16))
# plt = plt.ticklabel_format(useOffset=False, style='plain')
# fig, ax1 = plt.subplots()
# fig, ax1 = plt.subplots(figsize=(20,16))
fig, ax1 = plt.subplots(figsize=(20,16))

ax2 = ax1.twinx()
ax1.plot(x_list,y_list, color="black")
ax1.plot(x_list_pid,y_list_pid, color="red", alpha=0.5)
ax1.ticklabel_format(useOffset=False, style='plain')

ax2.plot(x_list_gc, y_list_gc, color="blue", alpha=0.5)
ax2.set_ylim([0,65]) 


ax2.ticklabel_format(useOffset=False, style='plain')



ax1.set_xlabel('Genome coordinates')
ax1.set_ylabel('# repeats', color="black")
ax2.set_ylabel('GC %', color='b')

# plt.show()
plt.savefig('Tm_pid_hits_plot_gc.png', dpi=300)








# PRINT PLOTTING DATA TO A CSV FILE
with open('plottting_data.csv', 'wb') as f:
    writer = csv.writer(f)
    # writer.writerows(izip(x_list, y_list_pid, y_list))
    writer.writerow(('Genome_coordinates', '70pid_hits', '10Tm_difference_hits'))
    writer.writerows(izip(x_list, y_list_pid, y_list))

    

    
    
    
    
    
    
    
############################################################
# Time to run the code: end timer
############################################################
t1 = time.time() 
total = t1-t0
total = ("{0:.2f}".format(round(total,2)))

print "Total time taken = ", total, "seconds"



