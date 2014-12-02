#!/usr/bin/env python

'''

Compiled by Felix Francis

Description:    Function to count the number of a pattern of a defined substrings in a given string.

Example:    Count(CGATATATCCATAG, ATA) is equal to 3 (not 2) since we should account for overlapping occurrences of pattern (ATA) in the given string CGATATATCCATAG.

'''

def pattern_count(string, pattern):
    count = 0
    for i in range(0,(len(string)-len(pattern))):
            if string[i:(len(pattern)+i)] == pattern:
                count = count+1
    print count


with open('../data/dataset_2_6.txt') as input_data:
    string_input,pattern = [line.strip() for line in input_data.readlines()]
    print pattern

pattern_count(string_input,pattern)


