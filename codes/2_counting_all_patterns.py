#!/usr/bin/env python

'''

Compiled by Felix Francis

Description:    Function to count the number of all frequent occuring k-mers in a given string.
                This approach is useful for smaller strings. When it comes to larger strings and greater number of k-mers, use this approach
                https://stepic.org/lesson/CS-The-Frequency-Array-2994/step/1

Example:    Count(ATAATATGAAGAAGAAGAAGCCATGATC, 3) gives AAG, GAA

'''

# Pattern count function
def pattern_count(string, pattern):
    count = 0
    for i in range(0,(len(string)-len(pattern))):
            if string[i:(len(pattern)+i)] == pattern:
                count = count+1
    return count

# Count the frequent k-mer
def frequent_words(string, k):
    frequent_patterns = set()
    count = 0
    for i in range(0,(len(string)-k)):
            pattern = string[i:(k+i)]
            # If the count of the k-mer is greater than the previous one the count is replaced by the new one, and the frequent_patterns set is replaced by the new frequent pattern
            if pattern_count(string,pattern) > count:
                count = pattern_count(string,pattern)
                frequent_patterns = set([pattern])
            # If the count of the k-mer is same as the previous one the new frequent pattern is added to the frequent_patterns set
            if pattern_count(string,pattern) == count:
                frequent_patterns.update([pattern])

    print frequent_patterns

##
#string1 = "ATAATATGAAGAAGAAGAAGCCATGATC"
#k_mer = 3
#frequent_words(string1, k_mer)
##

##
with open('../data/dataset_2_9.txt') as input_data:
    string_input,k_mer = [line.strip() for line in input_data.readlines()]
    #print type(k_mer)

frequent_words(string_input, int(k_mer))
##
