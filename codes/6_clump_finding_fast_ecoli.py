#!/usr/bin/env python

'''

Compiled by Felix Francis

Description:    A k-mer is a "clump" if it appears many times within a short interval of the genome.
                More formally, given integers L and t, a k-mer pattern forms an (L, t)-clump inside a (larger) string Genome,
                if there is an interval of Genome of length L in which this k-mer appears at least t times.
                
                Faster method using dictionary.

Input file:     A string Genome, and integers k, L, and t.

Example:    pattern_coord(CGATATATCCATAG, ATA) is 2, 4, 10 (keep in mind that the string and the out put coordinates starts from 0)
'''
# Check time to run the code
import time
t0 = time.time()


##Code start

dna = "CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA"
k   = 5 
L   = 50
t   = 4

#with open('../data/E-coli.txt') as input_data:
#    seq= [line.strip() for line in input_data.readlines()]
#    dna = seq[0]
#    k = 9
#    L = 500
#    t = 3
    

def CheckClumpLength(indicies, t, L):
	'''Checks that a given set of t k-mers falls within a clump of size L.'''
	for i in  xrange(len(indicies)-t+1):
		if indicies[t+i-1] - indicies[i-1] <= L:
			return True
	return False
	
#with open('../data/dataset_4_4.txt') as input_data:
#	dna, [k, L, t] = [line.strip() if index == 0 else map(int, line.strip().split()) for index, line in enumerate(input_data.readlines())]

# Find all k-mers, count their appearances, and store thier indicies. 
kmer_dict = dict()
for i in xrange(len(dna)-k):
	if dna[i:i+k] in kmer_dict:
		kmer_dict[dna[i:i+k]][0] += 1
		kmer_dict[dna[i:i+k]][1].append(i)
	else:
		kmer_dict[dna[i:i+k]] = [1, [i]]

# The candidate k-mers that appear at least t times, along with the indicies where they appear.
kmer_candidates = [ [kmer[0],kmer[1][1]] for kmer in kmer_dict.items() if kmer[1][0] >= t]

# Check that at least t candidate k-mers fall within a clump of size L.
kmer_clumps = []
for candidate in kmer_candidates:
	if CheckClumpLength(candidate[1], t, L):
		kmer_clumps.append(candidate[0])



# Print and save the solution.
print ' '.join(kmer_clumps)
print len(kmer_clumps)  
    
    
'''    
#print len(seq)
# create a dictionary of all k_mers, their frequencies and start coordinates
k_mer = dict()

for i in xrange(len(seq)-k+1):
    if seq[i:i+k] in k_mer:
        k_mer[seq[i:i+k]][0] +=1
        k_mer[seq[i:i+k]][1].append(i)
    else:
        k_mer[seq[i:i+k]] = [1,[i]]
#print k_mer

# Reduce the size of dictionary by selecting only the keys if their corresponding frequencies >= (t)
k_mer_candidates =[]
for k_mer in k_mer.items():
    if k_mer[1][0] >= t:
        k_mer_candidates.append([k_mer[0],[k_mer[1][0],k_mer[1][1]]])

# Checks whether a given set of t k-mers falls within a clump of size L
k_mer_clumps = set()
for candidate in (k_mer_candidates):
    for i in xrange(len(candidate [1][1])-t+1):
        if (candidate [1][1][t+i-1] - candidate [1][1][i]) <= L:
            k_mer_clumps.update([candidate [0]]) 
#print ' '.join(k_mer_clumps)
print len(k_mer_clumps)

##Code end
'''


t1 = time.time()

total = t1-t0
print total

