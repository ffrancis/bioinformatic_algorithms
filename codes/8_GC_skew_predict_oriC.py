#!/usr/bin/env python

'''

Compiled by Felix Francis

Description:    We can compute Skewi+1(Genome) from Skewi(Genome) according to the nucleotide in position i of Genome.
		If this nucleotide is G, then Skewi+1(Genome) = Skewi(Genome) + 1;
		if this nucleotide is C, then Skewi+1(Genome)= Skewi(Genome) - 1; otherwise, Skewi+1(Genome) = Skewi(Genome).
		
		Then find the index on the sequence corresponding to the smallest skew value, which would be the oriC position.

Input file:     A string Genome.

Example:    Input_string = "TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"
'''

from __future__ import print_function
# Check time to run the code
import time
t0 = time.time()

#Input file
with open ('../data/minimum_skew_data.txt') as input_data:
	dna = [line.strip() for line in input_data.readlines()]
	#dna = data[0]
	print (dna)


'''
dna = "TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"
'''

score  = [0]
def gc_skew_counter(sequence):
	sequence = sequence.upper()
	for i in sequence:
		if i == "G":
			score.extend([score[-1] +1])
		if i == "C":
			score.extend([score[-1] -1])
		if i != "G" and i != "C":
			score.extend([score[-1]])
	#print (*score, sep= ' ')
	#print (score)
	return score
		
gc_skew_counter(dna)


min_value, min_value_index = min((value, index) for (index, value) in enumerate(score))

smallest_indexes = []
for index, value in enumerate(score):
	if value == min_value:
		smallest_indexes.append(index)

print (*smallest_indexes, sep= ' ')

t1 = time.time()

total = t1-t0
print (total)

