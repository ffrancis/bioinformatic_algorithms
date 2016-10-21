### repeat quantification ofa genomic locus, for every non overlapping 25bp kmers.
### input files are the TA PSE module outputs (for exact match and no exact match options)
### quantifies the repeat count based on thermoalignment Tm and also based on number of mismatches


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
# from parameters import *
from itertools import repeat
import os.path


primer_size = 25

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
    allowed_chars    =    set('ACGT')
    if set(sequence.upper())    <=    allowed_chars:
        gc_count =  sequence.count('G')+sequence.count('C')
        sequence_len = len(sequence)
        gc_cont    = (float(gc_count)/sequence_len)*100
        gc_cont = round(gc_cont, 2)  
    else:
        gc_cont = np.nan
    return gc_cont  
     

# ADD ZEROS TO MISSING COORDINATES (USES start_pos, stop_pos FROM PARAMETER.PY)
def fill_in_missing_coords(x_list, y_list):
    for i in xrange(start_pos, stop_pos+1):
        if i in x_list:
            pass
        else:
            x_list.append(i)
            # y_list.append(0)
            y_list.append('-')

            
            # REMOVE UNDESIRED COORD INFO:
def remove_undesired_coords(x_list, y_list):
    desired_coords = range(start_pos, stop_pos, primer_size)
    for i in xrange(len(x_list)):
        if x_list[i] not in desired_coords:
            y_list[i] = '-'

            
            
### INPUT FILES AND PATHS

# file_path = './'
# directory_name = "qMDR_7_128386997_50394_repeats_viz_broad_param"         ### TO BE RUN
# directory_name = "qNLB_1_25376615_22184_repeats_viz_broad_param"          ### TO BE RUN
# directory_name = "qNLB_1_187278617_197947_viz_broad_param"                ### TO BE RUN
directory_name = "qSLB_2_37556845_13001_viz_broad_param"
# directory_name = "qSLB_3_219917184_72001_viz_broad_param"
# directory_name = "qSLB_6_7002788_135001_viz_broad_param"
# directory_name = "qSLB_9_16257370_95436_viz_broad_param"


from qSLB_2_37556845_13001_viz_broad_param.time_stamp import Time_stamp
from qSLB_2_37556845_13001_viz_broad_param.parameters import *

# from qMDR_7_128386997_50394_repeats_viz_broad_param.time_stamp import Time_stamp
# from qMDR_7_128386997_50394_repeats_viz_broad_param.parameters import *


# from qNLB_1_187278617_197947_viz_broad_param.time_stamp import Time_stamp
# from qNLB_1_187278617_197947_viz_broad_param.parameters import *


file_path_no_em = "/home/ffrancis/TA_codes_restored/LRPCR_050216/LRPCR_whole_loci/repeats_viz_100416/" + directory_name
file_path_em = file_path_no_em + "_em_picking"










pid_threshold = 70

f_primer_alignments_no_em = file_path_no_em + "/"+ Time_stamp+ "/"+ Time_stamp +"_HSE_out1_1.csv"
primer_coord_files_no_em  = file_path_no_em + "/"+ Time_stamp+ "/"+ Time_stamp +"_HSE_out3_1.csv"


f_primer_alignments_em = file_path_em + "/"+ Time_stamp+ "/"+ Time_stamp +"_HSE_out1_1.csv"
primer_coord_files_em  = file_path_em + "/"+ Time_stamp+ "/"+ Time_stamp +"_HSE_out3_1.csv"






### CODE



df_f_primers_em = pd.read_csv(f_primer_alignments_em)

df_f_primers_no_em = pd.read_csv( f_primer_alignments_no_em)


primer_names_no_em = pd.read_csv(primer_coord_files_no_em)
primer_names_em = pd.read_csv(primer_coord_files_em)






# MERGE EM AND NON-EM DATAFRAMES

df_f_primers_frames = [df_f_primers_no_em, df_f_primers_em]
primer_names_frames = [primer_names_no_em, primer_names_em]

df_f_primers = pd.concat(df_f_primers_frames)
primer_names = pd.concat(primer_names_frames)






# CREATE DICTIONARY OF PRIMER SEQUENCE AND PRIMER COORDINATE
primer_coord_dict = dict(zip(primer_names.Primer, primer_names.Genome_start))




# IF PRIMER LENGTH NOT EQUAL TO 25bp, REPLACE THAT PRIMER SEQUENCE WITH NaN
df_f_primers['Primer'] = np.where( df_f_primers['Primer'].str.len() == 25, df_f_primers['Primer'], np.nan)
# REMOVE ROWS WITH NaN AS PRIMER SEQUENCE (PRIMER LENGTHS != 25)
df_f_primers = df_f_primers.dropna(subset=['Primer']) 


'''
# IF PRIMER HAS NNN, REPLACE THAT PRIMER SEQUENCE WITH NaN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# s4.str.contains('A', na=False)
df_f_primers['Primer'] = np.where( df_f_primers['Primer'].str.contains('N'),  np.nan, df_f_primers['Primer'])
# REMOVE ROWS WITH NaN AS PRIMER SEQUENCE (PRIMER LENGTHS != 25)
df_f_primers = df_f_primers.dropna(subset=['Primer'])
'''



### FOR Tm BASED EVALUATION
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
 

        
fill_in_missing_coords(x_list, y_list_unsorted)
    
y_list = [x for (y,x) in sorted(zip(x_list, y_list_unsorted))]
x_list.sort()

remove_undesired_coords(x_list, y_list)


for i in xrange(0, len(y_list)-(primer_size), 25):
    if y_list[i] == '-':
        y_list[i] = 0
    # i = i +25



# for i in xrange(len(y_list)-25):
    # if y_list[i] == '-':
        # y_list[i] = 0
    # i = i +25





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

remove_undesired_coords(x_list_pid, y_list_pid)



for i in xrange(0, len(y_list_pid)-(primer_size), 25):
    if y_list_pid[i] == '-':
        y_list_pid[i] = 0
    # i = i +25





# GC CONTENT FROM WHOLE LOCUS

### Locus file name
locus       =    file_path_no_em +  "/" + Time_stamp +  "/" + Time_stamp + "_" + str(chr_no) + "_" + str(start_pos) + "_" + str((stop_pos - start_pos)+1) +"_VariantMasked.fasta"

### Process input loci sequence file
## If input file is .fasta file:
with open (locus) as sequence_data:
    line = sequence_data.read()
    lines = line.split("\n")
    input_seq = lines[1]
    


    
locus_gc_list = []
for i in xrange(len(input_seq)-primer_size+1):
    primer    = input_seq[i:i+primer_size].upper()
    primer_gc = gc_content(primer)
    locus_gc_list.append(primer_gc)
    
locus_gc_list.extend(repeat(np.nan, primer_size-1))
    
    
y_list_gc = locus_gc_list







'''
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
'''



# PRINT PLOTTING DATA TO A CSV FILE
with open(directory_name + '_repeats25bp_window101916.csv', 'wb') as f:
    writer = csv.writer(f)
    # writer.writerows(izip(x_list, y_list_pid, y_list))
    writer.writerow(('Genome_coordinates', '70pid_hits', '10Tm_difference_hits', 'GC'))
    writer.writerows(izip(x_list, y_list_pid, y_list, y_list_gc))

############################################################
# Time to run the code: end timer
############################################################
t1 = time.time() 
total = t1-t0
total = ("{0:.2f}".format(round(total,2)))

print "Total time taken = ", total, "seconds"



