#!/usr/bin/env python

'''

Compiled by Felix Francis (felixfrancier@gmail.com)

Description:    Predict DnaA boxes by identifying frequent k-mers, possibly with mismatches.
		Given strings Text and Pattern as well as an integer d, we define Countd(Text, Pattern) as the total number of
		occurrences of Pattern in Text with at most d mismatches. 
		
		Computing Countd(Text, Pattern) simply requires us to compute the Hamming distance between Pattern and every k-mer substring of Text


Example:    	For example, Count1(AACAAGCTGATAAACATTTAAAGAG, AAAAA) = 4
		because AAAAA appears four times in this string with at most one mismatch: AACAA, ATAAA, AAACA, and AAAGA.
		(Note that two of these occurrences overlap)
'''


from __future__ import print_function

#Time to run the code: start timer
import time
t0 = time.time()

'''
#Input file
with open ('../data/dataset_9_6.txt') as input_data:
	lines = input_data.readlines()
	data = [line.strip() for line in lines]
	#print (data [2])
	sequence = data[0]
	query = data[1]
	max_mismatch = int(data[2])	

'''

'''
#Test data

query = "GAGG"
sequence = "TTTAGAGCCTTCAGAGG"
max_mismatch = 2
'''


##Test2

query = "ATGT"
sequence = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
max_mismatch = 1
# Hamming distance function

def HammingDistance(p, q):
	score = 0
	if len(p) == len(q):
		for i in xrange (len(p)):
			if (p[i]) != q [i]:
				score = score+1
	return (score)


# Approximate pattern matching function

def approx_pattern_match(query, sequence, max_mismatch):
	index= []
	for i in xrange(len(sequence)-len(query)+1):
		pattern = sequence [i:i+len(query)]
		if HammingDistance(query,pattern) <= max_mismatch:
			index.append(i)
	#return ' '.join(map(str, index))
	return (index)
			
	

out_put_index = approx_pattern_match(query, sequence, max_mismatch)

print (len(out_put_index))


##Write the index to output.txt file
#with open('output.txt', 'w') as output_data:
#		output_data.write(out_put_index)
#



#Time to run the code: end timer
t1 = time.time()

total = t1-t0
print (total)

