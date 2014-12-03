#!/usr/bin/env python

'''

Compiled by Felix Francis

Question: Python script to identify the reverse compliment of a user input DNA sequence

Description:    Complementary strand of DNA contains the corresponding nucleotide bases as given in the dictionary (dict) in the script.
                But since DNA complementary strand runs antiparallel, we need to get its reverse (reverse compliment) to maintain the 5' - 3 ' orientation.
                
                For example:
                            Original Sequence  5'ATGCAGGGGAAACATGATTCAGGAC 3'
                            Complement         3'TACGTCCCCTTTGTACTAAGTCCTG 5'
                            Reverse Complement 5'GTCCTGAATCATGTTTCCCCTGCAT 3'

'''
dict = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}


'''
For user input sequence and validation of the DNA sequence:

import sys                                                                      # function to exit a script under if else conditions
sequence = raw_input("Please enter your sequence for getting it's reverse compliment: ")
sequence = sequence.upper()
print '\n',"The sequence you have entered is: ", '\t',  sequence
if sequence.count('A') + sequence.count('T') + sequence.count('G') + sequence.count('C') != len(sequence):
    print '\n',"Input DNA sequence should only contain the bases A,T,G or C: ",'\n','\n',
    sys.exit("Error message: Check your input sequence!")

sequence              = "ATGCAGGGGAAACATGATTCAGGAC"

'''


def rev_complement(seq):
    compliment_seq     = ''.join(dict.get(i) for i in seq)
    rev_compliment_seq = compliment_seq[::-1]
    print rev_compliment_seq


with open('../data/dataset_3_2.txt') as input_data:
    sequence = [line.strip() for line in input_data.readlines()]
    sequence = sequence[0]
    
rev_complement(sequence)

