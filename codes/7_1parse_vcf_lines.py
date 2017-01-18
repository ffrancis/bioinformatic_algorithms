### vcf_conversion: converts a vcf file to ThermoAlign input format
### Version 1.0.0: 06/28/2016
### Authors: Felix Francis (felixfrancier@gmail.com); Randall J. Wisser (rjw@udel.edu) 

### Requirements
### vcf files should follow the format described here: https://samtools.github.io/hts-specs/VCFv4.2.pdf
### all vcf input files must be named according to the corresponding fasta input files, e.g. chr1.vcf, chr2.vcf, ... 



### Import functions
import datetime
import os
import pandas as pd
import numpy as np

############################################################
### Time to run the code: start timer
############################################################
import time
t0 = time.time()


############################################################
### Functions
############################################################



input_path   =   "/mnt/data27/ffrancis/HapMap3/all_lines/"


select_line     = "Tx303"

### actual data
input_file   =   "c1_hmp31_q30.vcf"
locus_start = 25375814
locus_stop = 25398735

### test data
# input_file   =   "c10_first1000.vcf"
# locus_start =  229011
# locus_stop = 229703





df = pd.read_csv(input_path + input_file, sep='\t', comment='#', skiprows=0, usecols=[1, 4, 720], header=None)
# print df
# df = df[(df.coordinate >= int(start_pos)) & (df.coordinate <= int(stop_pos))]
df = df[(df[1] >= int(locus_start)) & (df[1] <= int(locus_stop))]


df.columns = ['coordinate', 'alternate_allele', "line"]
df.to_csv('parsed_snps' + input_file, sep='\t', encoding='utf-8', index=False)

