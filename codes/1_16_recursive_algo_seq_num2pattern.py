#!/usr/bin/env python

'''
Compiled by Felix Francis (felixfrancier@gmail.com)

Description:    in 1_15, we used the equation given below:

		PatternToNumber(AGT) = 4*PatternToNumber(AG) + SymbolToNumber(T) 

		Now, if we need to do the reverse, i.e., compute the inverse function NumberToPattern(index, k), based on the equation
		given above, when we divide index = PatternToNumber(Pattern) by 4, the remainder will be equal to SymbolToNumber(symbol),
		and the quotient will be equal to PatternToNumber(Prefix(Pattern)). Thus, we can use this fact to peel away symbols at
		the end of Pattern one at a time as illulstrated in the link: http://www.ucarecdn.com/229ad845-ca9d-450f-a03c-cfb30d338908/
		
Example:	number2pattern(11,3) gives AGT

'''
############################################################
#Time to run the code: start timer
############################################################
import time
t0 = time.time()


############################################################
# Function to get the base(last symbol) based on the last number
# eg. bases = "ACGT" correspond to 0123; so 1 = "C"
############################################################
def number2symbol(last_number):
	bases = "ACGT"
	if last_number >=0 and last_number <=3:
		return bases[last_number]

#print (number2symbol(1))


############################################################
# Function to get the pattern based on the index number and "k" pattern size
############################################################
def number2pattern(number,k):
	if isinstance(number, int) and isinstance(k, int):
		if k == 1:
			return number2symbol(number)
		quotient = number/4
		remainder = number%4
		pattern = number2pattern(quotient,k-1)
		symbol = number2symbol(remainder)
		return (str(pattern)+ str(symbol))	
	else:
		return "Make sure your index number and pattern size are integers!!"

print(number2pattern(6426,7))





############################################################
#Time to run the code: end timer
############################################################
t1 = time.time()

total = t1-t0
print '\n', "Time to run code = ", total, " seconds", '\n'




















