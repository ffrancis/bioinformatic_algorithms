#!/usr/bin/env python

'''
Compiled by Felix Francis (felixfrancier@gmail.com)

Description:    1. Compute the effective cnv space (ecs) in a sequence.
                2. Within the effective cnv space (ecs), identify the variants (secondary) and their coordinates within the ecs.
                3. Similar variants (secondary) in different sequences in the MSA (if any variant is not present in the same position) are identified as the same,
                   since they should have similar biological impacts.

Example:
            Input FASTA sequences:
		>Sequence1:
		ATGCTXTGTATGCCGTTAAG
		>Sequence2:
		ATGCYYTGTATGCCGTTAAG 
		>Sequence3:
		ATGCTXTGTATGCCGTTAAG 
		>Sequence4:
		ATGCTXTGTATGTXTGCCGTTAAG
		>Sequence5:
		ATGCTATGTXTGTXTGCCGTTAAG

            ClustalW out put:
		Sequence1_      ATGC TXTG TATG ---- CCGT TAAG 
		Sequence3_      ATGC TXTG TATG ---- CCGT TAAG 
		Sequence4_      ATGC TXTG TATG TXTG CCGT TAAG 
		Sequence5_      ATGC TATG TXTG TXTG CCGT TAAG 
		Sequence2_      ATGC YYTG TATG ---- CCGT TAAG 


            True alignment should have been:
		Sequence1_      ATGC TXTG TATG ---- ---- CCGT TAAG 
		Sequence3_      ATGC TXTG TATG ---- ---- CCGT TAAG 
		Sequence4_      ATGC TXTG TATG TXTG ---- CCGT TAAG 
		Sequence5_      ATGC ---- TATG TXTG TXTG CCGT TAAG 
		Sequence2_      ATGC YYTG TATG ---- ---- CCGT TAAG 

                
            Because of the preferred left directed alignment of cnv's, variants in the segments constituting cnv's may be placed in different positions.
            This leads to wrong annotation of the variants (TXTG with X indicating a variation from TATG) in the segments constituting cnv's.
            When in reality they may have the same biological effects.

            Considering this, all the (TXTG) in different positions within the effective cnv polymorphism space (ecps) will be identified as the same unless there are
            TXTG in the same position in different sequences or if there are more than one

Notes:
            Consider the largest of the cnv's if there are smaller repats/cnv's within the largest ones.
	    In future, include cnvs of different lengths.
	    
	    same index as query, pattern index1 + key, pattern index2 + key,......
	    
	    What if YY in Sequence2_ are deletions. Then "-" should not be included while using this algorithm
'''
############################################################
#Time to run the code: start timer
############################################################
import time
t0 = time.time()


# Import all the required functions
import numpy


############################################################
## Assign parameters
############################################################
# min_cnv_size is the smallest size of cnv's to be checked
min_cnv_size = 3
# percentage identity between sequences to be counted as cnv's
p_id = 80  # Use lower p_id values < 70  for cnv's constituted by smaller sequences (>9bp ) for larger sequences that make up cnv's use higher p_id's


############################################################
#Input data
############################################################
seq_name = []
seq = []
#with open ('data/clustalw2-I20150106-042737-0540-50781400-oy.clustalw') as input_data:
with open ('data/test.clustalw') as input_data:	
	lines1 = input_data.readlines()	
	for line in lines1:
		line2=line.split('\n')
		cols = line2[0].split(' ')
		if cols[0] != '' and cols[0] != "CLUSTAL":
			seq_name.append(cols[0] [:-1])
			seq.append(cols[6])

#print seq
############################################################
# Search cnv's
############################################################

# Function to compare query with target and return the number of mismatches
def seq_compare(p, q):
	score = 0.0
	if len(p) == len(q):
		for i in xrange (len(p)):
			if (p[i]) != q [i]:
				score = score+1
	return (score)
	#print score


# Pattern matching function for identifying cnv's which are one after another

###############################################################################
##for s in xrange(len(seq)):
##	sequence = seq[s].replace("-", "")
##	print sequence
##	# Get the length of aligned sequences (ith one)
##	seq_len = len(sequence) # in the code below, change seq [3] to sequence
###############################################################################

# Get the length of aligned sequences (ith one)
seq_len = len(seq[3])  					# replace seq [3] with seq [s]

# List of k_mers sizes greater than min_cnv_size upto the length of the aligned sequence.
k_mer_sizes = xrange(min_cnv_size,(seq_len/2)+1)
array = numpy.array(k_mer_sizes)
ext_space = seq_len - array

# Dictionary of k_mer_sizes as keys and corresponding distances (as list) to be covered while searching for cnv's
k_mer_sizes_dict = dict((key,0) for key in k_mer_sizes)
for key in k_mer_sizes_dict.keys():
	iter_coord_list = range(0, seq_len - key + 1)
	k_mer_sizes_dict[key] = iter_coord_list




index1 = []
for key, value in k_mer_sizes_dict.iteritems():
	query = []
	patterns = []
	max_mismatch = []
	#print key,(value)
	for i in value[:-key]:
		query.append(seq[3] [i:i+(key)])					# replace seq [3] with seq [s] / to sequence (no "-") if iterating over all sequence
		patterns.append(seq[3] [(i+key):((i+key)+(key))])			# replace seq [3] with seq [s] / to sequence (no "-") if iterating over all sequence
		max_mismatch.append(round( ((100.0-p_id)/100.0)*(len(seq[3] [i:i+(key)]))) )	# replace seq [3] with seq [s] / to sequence (no "-") if iterating over all sequence
		#max_mismatch.append(( ((100.0-p_id)/100.0)*(len(seq[3] [i:i+(key)]))) )	# replace seq [3] with seq [s] / to sequence (no "-") if iterating over all sequence

	#print "Query = ", (query)
	#print "Patterns = ", (patterns)
	#print "max_mismatch = ", max_mismatch 
	#print "# Query = ", len(query)
	#print "# Patterns = ", len(patterns)

	
	
	cnv_dict = dict()
	#index = []
	for j in xrange(len(query)):
		#print j
		if seq_compare(query[j],patterns [j]) <= max_mismatch [j]:		
			#print "True"
			#index.append(j)				# counting starts from 0!!!
			#index.append(j+key)
			cnv_dict[query[j]] = [j, [j+key]]
			#for k in xrange(len(patterns)):	
			#	if k+key < len(patterns):
			#		if seq_compare(query[j],patterns [k+key]) < max_mismatch [j]:		
			#			print "True"
			#			index.append(k+key)
						

	#print "cnv index = ", index
	print cnv_dict


#Time to run the code: end timer
t1 = time.time()

total = t1-t0
print '\n', "Time to run code = ", total, " seconds", '\n'




















