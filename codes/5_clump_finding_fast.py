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

#seq = "ATAATAAATGATAATAAATG"
#k   = 3 
#L   = 6
#t   = 2

with open('../data/dataset_4_4.txt') as input_data:
    seq,k_L_t= [line.strip() for line in input_data.readlines()]
    k_L_t = (k_L_t.split())
    k = int(k_L_t[0])
    L = int(k_L_t[1])
    t = int(k_L_t[2])

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
print ' '.join(k_mer_clumps)

##Code end



t1 = time.time()

total = t1-t0
print total

