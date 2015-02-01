#!/usr/bin/env python

'''
Compiled by Felix Francis (felixfrancier@gmail.com)

Description:    If we remove the final symbol from all lexicographically ordered k-mers, the resulting list is still ordered
		lexicographically (think about removing the final letter from every word in a dictionary). In the case of DNA strings,
		every (k-1)-mer in the resulting list is repeated four times. Thus, the number of 3-mers occurring before AGT is
		equal to four times the number of 2-mers occurring before AG plus the number of 1-mers occurring before T. Therefore,

		PatternToNumber(AGT) = 4*PatternToNumber(AG) + SymbolToNumber(T) = 8 + 3 = 11

		where SymbolToNumber(symbol) is the function transforming symbols A, C, G, and T into the respective integers 0, 1, 2, and 3.

Example:

'''
############################################################
#Time to run the code: start timer
############################################################
import time
t0 = time.time()

############################################################
# Function to get the number of bases beofore a given base, when the bases are ordered lexicographically
# (eg. symbol2number("T") == 2 since in a and C in"ACTG" are before "T")
############################################################
def symbol2number(symbol):
	symbol= symbol.upper()
	bases = "ACGT"
	if set(symbol) <= set(bases):
		return bases.index(symbol)
		#print bases.index(symbol)

#symbol = "g"
#symbol2number(symbol)

############################################################
# Recursive function to count the number of k-mers occurring before a given pattern
############################################################
def pattern2number(pattern):
	pattern = pattern.upper()
	bases = set("ATGC")
	if set(pattern) > bases:
		print "Wrong bases in input sequence pattern"
	if len(pattern) == 0:
		return 0
	else:
		symbol = pattern[-1:]
		pattern = pattern[:-1]
		return (4*(pattern2number(pattern))) + symbol2number(symbol)

	
pattern = "AGT"
print(pattern2number(pattern))


############################################################
#Time to run the code: end timer
############################################################
t1 = time.time()

total = t1-t0
print '\n', "Time to run code = ", total, " seconds", '\n'




















