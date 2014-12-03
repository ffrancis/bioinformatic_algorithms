#!/usr/bin/env python

'''

Compiled by Felix Francis (felixfrancier@gmail.com)

Description:    Position i in k-mers pi .... pk and qi .... qk is a mismatch if pi != qi.
		For example, CGAAT and CGGAC have two mismatches.
		The number of mismatches between strings p and q is called the Hamming distance between these strings
		and is denoted HammingDistance(p, q).


Example:    	Input_string p = "GGGCCGTTGGT"
		Input_string q = "GGACCGTTGAC"
'''


from __future__ import print_function

#Time to run the code: start timer
import time
t0 = time.time()


#Input file
with open ('../data/dataset_9_4.txt') as input_data:
	lines = input_data.readlines()
	data = [line.strip() for line in lines]
	#print (data [2])
	query = data[0]
	sequence = data[1]
	max_mismatch = int(data[2])	
	


#Test data
'''
query = "ATTCTGGA"
sequence = "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT"
max_mismatch = 3
'''


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
	#print type(index)
	return ' '.join(map(str, index))
			
	

out_put_index = approx_pattern_match(query, sequence, max_mismatch)

#Write the index to output.txt file

with open('output.txt', 'w') as output_data:
		output_data.write(out_put_index)




#Time to run the code: end timer
t1 = time.time()

total = t1-t0
print (total)

