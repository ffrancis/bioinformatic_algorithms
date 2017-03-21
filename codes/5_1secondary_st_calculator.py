
############################################################
#Time to run the code: start timer
############################################################
import time
t0 = time.time()


'''

# INPUT DATA FORMAT:
11.4_01_f	GTTAGTCAGACGATGCGTCATGCGGCTGTGAAGATGTCATGGACA	11.4_01_r	GTTAGATGACGCATCGTCTGATGGAACAAGAGCTGGCCTGTTCGA
11.4_02_f	GTTAGCTATACATGACTCTGCGCGGCTGTGAAGATGTCATGGACA	11.4_02_r	GTTAGGCAGAGTCATGTATAGTGGAACAAGAGCTGGCCTGTTCGA
11.4_03_f	GTTAGTACTAGAGTAGCACTCGCGGCTGTGAAGATGTCATGGACA	11.4_03_r	GTTAGGAGTGCTACTCTAGTATGGAACAAGAGCTGGCCTGTTCGA

'''




############################################################
#### IMPORT FUNCTIONS
############################################################
import math
import re
from santalucia_tm_model091615 import NN_Tm                
from santalucia_tm_model091615 import complement        
from santalucia_tm_model091615 import mM_monovalent        
import csv
import primer3
# import subprocess as sp                                   
# from operator import itemgetter
# import shutil
# import os.path
# import glob
from parameters import *

### FUNCTION TO CALCULATE HAIRPIN TM
def hairpin_Tm(primer_sequence, mv_cation=0,primer_conc=0): 
    Tm_hairpin =  (primer3.calcHairpin(primer_sequence,mv_conc=mv_cation, dv_conc=0, dntp_conc=0, dna_conc=primer_conc, temp_c=37, max_loop=30)).tm
    return ("{0:.2f}".format(round(Tm_hairpin,2)))    

### FUNCTION TO CALCULATE HOMODIMER TM
def homodimer_Tm(primer_sequence, mv_cation=0,primer_conc=0): 
    Tm_homodimer = (primer3.calcHomodimer(primer_sequence,mv_conc=mv_cation, dv_conc=0, dntp_conc=0, dna_conc=primer_conc, temp_c=37, max_loop=30)).tm
    return ("{0:.2f}".format(round(Tm_homodimer,2)))        
    
### function to check hetero dimer dg
def heterodimer_dg(seq1, seq2, mv_cation=0,primer_conc=0):
        dg =  (primer3.calcHeterodimer(seq1, seq2,mv_conc=mv_cation, dv_conc=0, dntp_conc=0, dna_conc=primer_conc, temp_c=60, max_loop=30)).tm
        return float(("{0:.2f}".format(round(dg,2))))

        
### Compute the monovalent cation equivalent
monovalent_cation_eq    =    mM_monovalent(Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs)


output_sec_st_info = open("secondary_structure_info.txt", "w")
output_sec_st_info.write('fPrimer_name'+'\t'+'fPrimer_seq'+'\t'+'fPrimer_Tm_hairpin'+'\t'+'fPrimer_homodimer'+'\t'+'rPrimer_name'+'\t'+'rPrimer_seq'+'\t'+'rPrimer_Tm_hairpin'+'\t'+'rPrimer_homodimer'+'\t'+'heterodimer_dg' +'\n')
input_file  =   "barcoded_primers.txt"
f = open(input_file)
csv_f = csv.reader(f, delimiter='\t')
for row in csv_f:
    f_primer, r_primer = row[1], row[3]
    f_primer_name, r_primer_name = row[0], row[2]
    print row[1], row[3]
    # print heterodimer_dg(f_primer, r_primer, monovalent_cation_eq, primer_conc)
    fTm_hairpin      =    str(hairpin_Tm(f_primer, monovalent_cation_eq, primer_conc))
    rTm_hairpin      =    str(hairpin_Tm(r_primer, monovalent_cation_eq, primer_conc))
    fTm_homodimer    =    str(homodimer_Tm(f_primer, monovalent_cation_eq, primer_conc))
    rTm_homodimer    =    str(homodimer_Tm(r_primer, monovalent_cation_eq, primer_conc))
    # heterodimer_dg  =    heterodimer_dg(f_primer, r_primer, monovalent_cation_eq, primer_conc)

    # print fTm_hairpin, rTm_hairpin, fTm_homodimer, rTm_homodimer, heterodimer_dg
    print fTm_hairpin, rTm_hairpin, fTm_homodimer, rTm_homodimer, heterodimer_dg(f_primer, r_primer, monovalent_cation_eq, primer_conc)
    output_sec_st_info.write(f_primer +'\t'+ f_primer_name +'\t'+fTm_hairpin+'\t'+fTm_homodimer+'\t'+r_primer +'\t'+ r_primer_name +'\t'+rTm_hairpin+'\t'+rTm_homodimer+'\t'+str(heterodimer_dg(f_primer, r_primer, monovalent_cation_eq, primer_conc)) +'\n')

output_sec_st_info.close() 
    
############################################################
#Time to run the code: end timer
############################################################
t1 = time.time()
total = t1-t0
print '\n', "Time to run code = ", total, " seconds", '\n'