#!/usr/bin/env python

'''

Compiled by Felix Francis

Description:    Position i in k-mers pi .... pk and qi .... qk is a mismatch if pi != qi.
		For example, CGAAT and CGGAC have two mismatches.
		The number of mismatches between strings p and q is called the Hamming distance between these strings
		and is denoted HammingDistance(p, q).


Example:    	Input_string p = "GGGCCGTTGGT"
		Input_string q = "GGACCGTTGAC"
'''


from __future__ import print_function
# Check time to run the code
import time
t0 = time.time()


#Input file
with open ('../data/dataset_9_3.txt') as input_data:
	lines = input_data.readlines()
	sequences = [line.strip() for line in lines]
	#print (sequences)
	p = sequences[0]
	q = sequences[1]
	
	
##Code start

'''
#Test data

p = "GGGCCGTTGGT"
q = "GGACCGTTGAC"
'''

# Hamming distance function

def HammingDistance(p, q):
	score = 0
	if len(p) == len(q):
		for i in xrange (len(p)):
			if (p[i]) != q [i]:
				score = score+1
	print (score)

HammingDistance(p, q)


t1 = time.time()

total = t1-t0
print (total)

