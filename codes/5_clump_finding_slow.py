#!/usr/bin/env python

'''

Compiled by Felix Francis

Description:    A k-mer is a "clump" if it appears many times within a short interval of the genome.
                More formally, given integers L and t, a k-mer pattern forms an (L, t)-clump inside a (larger) string Genome,
                if there is an interval of Genome of length L in which this k-mer appears at least t times.
                
                Slower method. Suitable for small sized sequences!!!

Input file:     A string Genome, and integers k, L, and t.

Example:    pattern_coord(CGATATATCCATAG, ATA) is 2, 4, 10 (keep in mind that the string and the out put coordinates starts from 0)
'''

import time

t0 = time.time()
#code_block




# Pattern count function
def pattern_count(string, pattern):
    count = 0
    for i in xrange(0,(len(string)-len(pattern))):
            if string[i:(len(pattern)+i)] == pattern:
                count = count+1
    return count

##
    
#string = "ATAATAAATGATAATAAATG"
#k   = 3
#L   = 10
#t   = 2
#
#clump_pattern_count(string,k,L,t) (take this below the function given below)


def clump_pattern_count(string,k,L,t):
    clump_pattern = set()
    for i in xrange(0,(len(string)-L+1)):
        sub_string = string[i:(i+L)]
        
        frequent_patterns = set()
        for i in xrange(0,(len(sub_string)-k)):
            pattern = sub_string[i:(k+i)]
            # If the count of the k-mer is greater than or same as t, the frequent pattern is added to the frequent_patterns set
            if pattern_count(sub_string,pattern) >= t:
                frequent_patterns.update([pattern])
    print frequent_patterns

###
with open('../data/dataset_4_4.txt') as input_data:
    string,k_L_t= [line.strip() for line in input_data.readlines()]
    k_L_t = (k_L_t.split())
    k = int(k_L_t[0])
    L = int(k_L_t[1])
    t = int(k_L_t[2])

##
clump_pattern_count(string,k,L,t)
###



t1 = time.time()

total = t1-t0
print total

