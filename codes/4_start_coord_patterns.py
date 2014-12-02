#!/usr/bin/env python

'''

Compiled by Felix Francis

Description:    Function to identify coordinates of defined substrings in a given string.

Example:    pattern_coord(CGATATATCCATAG, ATA) is 2, 4, 10 (keep in mind that the string and the out put coordinates starts from 0)
'''

def pattern_coord(string, pattern):
    start_coord = []
    for i in range(0,(len(string)-len(pattern))):
            if string[i:(len(pattern)+i)] == pattern:
                start_coord.extend([(i)])
    print start_coord

#string_input = "CGATATATCCATAG"
#pattern = "ATA"
pattern = "CTTGATCAT"

#with open('../data/dataset_3_5.txt') as input_data:
#    pattern, string_input = [line.strip() for line in input_data.readlines()]
#    print pattern
  
with open('../data/Vibrio_cholerae.txt') as input_data:
    string_input = [line.strip() for line in input_data.readlines()]
    string_input = string_input[0]  
    
pattern_coord(string_input,pattern)
