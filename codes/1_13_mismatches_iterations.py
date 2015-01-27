#!/usr/bin/env python

'''

Compiled by Felix Francis (felixfrancier@gmail.com)

Description:    A most frequent k-mer with up to d mismatches in Text is simply a string Pattern maximizing Countd(Text, Pattern)
		among all k-mers. Note that Pattern does not need to actually appear as a substring of Text.
		
Example:	AAAAA is the most frequent 5-mer with 1 mismatch in AACAAGCTGATAAACATTTAAAGAG,
		even though it does not appear exactly in this string.

'''


from __future__ import print_function

#Time to run the code: start timer
import time
t0 = time.time()


##Input file
#with open ('../data/dataset_9_7.txt') as input_data:
#	lines = input_data.readlines()
#	data = [line.strip() for line in lines]
#	#print (data)
#	sequence = data[0]
#	sub_data =str(data[1])
#	sub_data = [int(s) for s in sub_data.split()]
#	k = int(sub_data[0])
#	max_mismatch = int(sub_data[1])
#	print (max_mismatch )
#



##
##Test data
#
sequence 	= "ACGTTGCATGTCGCATGATGCATGAGAGCT"
k		= 4					#length of k-mers (query)
max_mismatch 	= 2




#sequence 	= "CACAGTAGGCGCCGGCACACACAGCCCCGGGCCCCGGGCCGCCCCGGGCCGGCGGCCGCCGGCGCCGGCACACCGGCACAGCCGTACCGGCACAGTAGTACCGGCCGGCCGGCACACCGGCACACCGGGTACACACCGGGGCGCACACACAGGCGGGCGCCGGGCCCCGGGCCGTACCGGGCCGCCGGCGGCCCACAGGCGCCGGCACAGTACCGGCACACACAGTAGCCCACACACAGGCGGGCGGTAGCCGGCGCACACACACACAGTAGGCGCACAGCCGCCCACACACACCGGCCGGCCGGCACAGGCGGGCGGGCGCACACACACCGGCACAGTAGTAGGCGGCCGGCGCACAGCC"
#k		= 7					#length of k-mers (query)
#max_mismatch 	= 2 


pattern = "AAGTC"

'''
AA

TA GA CA
AT AG AC
TT TG TC
GT GG GC
CT CG CC
GT GG GC
'''

bases = ['A','T','G','C']

variants = {}

'''
Generate max_mismatch variants of k-mers

from copy import deepcopy

ch = ['a','t','c','g']
d = {}
t = ['a', 't', 't']

def compute_replacements(s):
    for i in xrange(len(s)):
        for c in ch:
            if c!=s[i]:
                new_t = deepcopy(s)
                new_t[i] = c
                if tuple(new_t) not in d:
                    d[tuple(new_t)] = 0

d[tuple(t)] = 0
compute_replacements(t)
new_d = deepcopy(d)
for r in new_d:
    compute_replacements(list(r))
for r in d:
    print list(r)

'''




#Time to run the code: end timer
t1 = time.time()

total = t1-t0
print (total)

