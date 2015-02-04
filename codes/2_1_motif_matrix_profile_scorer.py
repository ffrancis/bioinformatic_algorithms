#!/usr/bin/env python

'''
Compiled by Felix Francis (felixfrancier@gmail.com)

Description:	1) Predict the approximate concensus motif (crude method) from a given set of sequence motifs
		2) Calculate the score of the given set of motifs

		Counts the number of occurrences of each nucleotide in each column of the motif matrix;
		the (i, j)-th element of Count(Motifs) stores the number of times that nucleotide i appears
		in column j of Motifs. We will further divide all of the elements in the count matrix by t,
		the number of rows in Motifs. This results in a profile matrix P = Profile(Motifs) for which Pi,j
		is the frequency of the i-th nucleotide in the j-th column of the motif matrix.
		Note that the elements of any column of the profile matrix sum to 1.
		
		Concensus motif from a given set of sequence motifs is identified based on the maximum profile score.
		
Example:	T   C   G   G   G   G   g   T   T   T   t   t           
		c   C   G   G   t   G   A   c   T   T   a   C
		a   C   G   G   G   G   A   T   T   T   t   C
		T   t   G   G   G   G   A   c   T   T   t   t
		a   a   G   G   G   G   A   c   T   T   C   C
		T   t   G   G   G   G   A   c   T   T   C   C
		T   C   G   G   G   G   A   T   T   c   a   t
		T   C   G   G   G   G   A   T   T   c   C   t
		T   a   G   G   G   G   A   a   c   T   a   C
		T   C   G   G   G   t   A   T   a   a   C   C
		
Score:		3 + 4 + 0 + 0 + 1 + 1 + 1 + 5 + 2 + 3 + 6 + 4 = 30

Count      A:   2   2   0   0   0   0   9   1   1   1   3   0          
	   C:   1   6   0   0   0   0   0   4   1   2   4   6  
	   G:   0   0  10  10   9   9   1   0   0   0   0   0  
	   T:   7   2   0   0   1   1   0   5   8   7   3   4

Profile:
	   A:  .2  .2   0   0   0   0  .9  .1  .1  .1  .3   0            
	   C:  .1  .6   0   0   0   0   0  .4  .1  .2  .4  .6  
	   G:   0   0   1   1  .9  .9  .1   0   0   0   0   0  
	   T:  .7  .2   0   0  .1  .1   0  .5  .8  .7  .3  .4

Consensus       T   C   G   G   G   G   A   T   T   T   C   C

'''
############################################################
#Time to run the code: start timer
############################################################
import time
t0 = time.time()

import pandas as pd



#Input file
with open ('../data/motif_matrix.txt') as input_data:
	seq = [line.strip() for line in input_data.readlines()]
	seq1 =[]
	for i in seq:
		j = (str(i).upper()).split()
		seq1.append(j)
		
#print len(seq1[0])		

df = pd.DataFrame(seq1)

#print (df[0]).str.contains(r'A').sum()

############################################################
#List count score function
############################################################
def list_count_char_score(list):
	lcDict= {}
	for i in list:
		if i in lcDict:
			lcDict[i] +=1
		else:
			lcDict[i]= 1
	#return lcDict
	import operator
	return (len(list)-max(lcDict.iteritems(), key=operator.itemgetter(1))[1])

#list= df[0].values.tolist()
#print(list_count_char_score(list))

############################################################
#List count max count base function
############################################################
def list_count_base(list):
	lcDict= {}
	for i in list:
		if i in lcDict:
			lcDict[i] +=1
		else:
			lcDict[i]= 1
	#return lcDict
	import operator
	return (max(lcDict.iteritems(), key=operator.itemgetter(1))[0])

#list= df[0].values.tolist()
#print(list_count_base(list))

############################################################
#Parse colums in Pandas data frame and calculate the total 
############################################################
score_list = []
concensus_motif =""
for i in (xrange(len(df.columns))):
	list= df[i].values.tolist()
	score_list.append(list_count_char_score(list))
	concensus_motif+=str(list_count_base(list))
final_score = sum(score_list)


print "Consensus motif = ", concensus_motif
print "Score of given set of motifs = ", final_score

############################################################
#Time to run the code: end timer
############################################################
t1 = time.time()

total = t1-t0
print '\n', "Time to run code = ", total, " seconds", '\n'




















