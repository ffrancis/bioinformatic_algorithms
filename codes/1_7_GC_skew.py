#!/usr/bin/env python

'''

Compiled by Felix Francis

Description:    We can compute Skewi+1(Genome) from Skewi(Genome) according to the nucleotide in position i of Genome.
		If this nucleotide is G, then Skewi+1(Genome) = Skewi(Genome) + 1;
		if this nucleotide is C, then Skewi+1(Genome)= Skewi(Genome) - 1; otherwise, Skewi+1(Genome) = Skewi(Genome).

Input file:     A string Genome.

Example:    Input_string = "GAGCCACCGCGATA"
'''

from __future__ import print_function
# Check time to run the code
import time
t0 = time.time()



##Code start

dna = "GAGCCACCGCGATA"

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
	print (*score, sep= ' ')
		
gc_skew_counter(dna)
	





t1 = time.time()

total = t1-t0
print (total)

