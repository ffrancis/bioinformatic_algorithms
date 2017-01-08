### Demultiplex and assemble error corrected PacBio reads. 
### Error corrected reads from each barcode (separate directory) and from each locus are pooled and assembled.
### Version 0.2.0: 12/19/2016
### Author: Felix Francis (felixfrancier@gmail.com)

############################################################
#Time to run the code: start timer
############################################################
import time
t0 = time.time()

############################################################
#### IMPORT FUNCTIONS
############################################################

import multiprocessing
from Bio import SeqIO
import pandas as pd
import datetime
import re
import os
from pydna import Assembly, Dseqrecord
import subprocess as sp  
import os.path
import time
import shutil


############################################################
#### PATHS AND INPUT FILES
############################################################


primer_pairs    = "primer_pairs_info.txt"
barcode_list    = "smrt1_barcodes_list.txt"
amos_path       = "/usr/local/amos/bin/"
consensus_output = "./output/"
trim_bp         = 21                           ### number of bases corresponding to padding + barcode that need to be trimmed from the amplicon consensus

### parameters
no_reads_threshold  = 100
cores               = 0        
nthread             = 3

############################################################
#### FUNCTIONS
############################################################

### pool reads from each locus
def pool_reads(primer_info_file, locus_chr_no, locus_start_pos, locus_stop_pos, locus_directory, required_barcode):
    df_primers = pd.read_csv(primer_info_file, sep='\t', skiprows=0, header=0)  ### read primer info
    reads = locus_directory +"reads_" + str(required_barcode)
    merged_fasta = open(reads +".fasta", "w")
    for index, row in df_primers.iterrows():
        f_primer_name, r_primer_name =  str(row['f_primer_name']), str(row['r_primer_name'])
        primer_chr_no, amplicon_start, amplicon_stop = int(f_primer_name.split("_")[1]), int(f_primer_name.split("_")[2]), int(r_primer_name.split("_")[2])
        if primer_chr_no == locus_chr_no and amplicon_start >= locus_start_pos and amplicon_stop <= locus_stop_pos:
            amplicon_start_stop_list = [amplicon_start, amplicon_stop]
            seq_path = "./amplicon" + "_" + str(locus_chr_no) + "_" +  str(amplicon_start) + "_" + str(amplicon_stop) + "/"
            fasta_sequences = SeqIO.parse(open(seq_path + "amplicon_analysis.fasta"), 'fasta') 
            for fasta in fasta_sequences:
                header, sequence = fasta.id.split('_'), str(fasta.seq)
                barcode = int(header[0][7:])    ### change to integer when acutal barcode data is available???                                                                                               
                no_reads = int(header[3][8:])
                id = seq_path.split("/")[1]+"_"+str(barcode)
                if barcode == required_barcode and no_reads >= no_reads_threshold:
                    merged_fasta.write(">" +str(barcode) + "_" +str(primer_chr_no) + "_" +str(amplicon_start)+ "_" +str(amplicon_stop) + "\n")
                    merged_fasta.write(sequence + "\n")
    merged_fasta.close()
    
### test pool_reads function
#pool_reads(primer_pairs,locus_chr_no, locus_start_pos, locus_stop_pos)



 


### convert a given fasta file containing sequence reads to afg format and carry out minimus assembly.
def minimus_assembly(barcode_directory, amos_path, consensus_output):
    ### convert fasta to afg format
    reads = consensus_output + barcode_name + "/" + barcode_name + "merged_reads"
    p = sp.Popen(["%stoAmos" %amos_path,"-s","%s" %reads+".fasta", "-o", "%s" %reads+"_assembly.afg"], stdout=sp.PIPE)
    out, err = p.communicate()
    ### assemble afg reads
    p = sp.Popen(["%sminimus" %amos_path, "%s" %reads+"_assembly.afg"], stdout=sp.PIPE)
    out, err = p.communicate()
    shutil.rmtree(reads+ "_assembly.bnk")
    os.remove(reads+ "_assembly.afg.runAmos.log")
    os.remove(reads+ "_assembly.afg") 
    os.remove(reads+ "_assembly.contig") 
    
    
    
### pick reads from each locus and assemble 
def pooled_locus_read_assembly(loci_info, primer_pairs, barcodes_list):
    ### create output parent directory
    out_parent_dir  = "./TESTAssembly_output_" + str((datetime.datetime.now().isoformat()).replace(":", "_").replace(".", "_")) +"/"
    os.makedirs(out_parent_dir)
    
    ### tab .... for each barcode combination
    
    
    ### read loci info
    df_loci = pd.read_csv(loci_info, sep='\t', skiprows=0, header=0)    ### read loci info, create directories for each locus
    for index, row in df_loci.iterrows():
        locus_name, locus_chr_no, locus_start_pos, locus_stop_pos = row['locus_name'], row['chr_no'], row['start_pos'], row['stop_pos']

        ### create directory for each locus (if it does not already exist)
        locus_directory = out_parent_dir +  str(locus_name)+ "/"
        if not os.path.exists(locus_directory):
            os.makedirs(locus_directory)
        ### pool_reads
        pool_reads(primer_pairs, locus_chr_no, locus_start_pos, locus_stop_pos, locus_directory, required_barcode)
        ### assemble reads from a locus
        minimus_assembly(locus_directory, required_barcode, amos_path)



############################################################
#### CODE
############################################################


### iterate barcode

### iterate primer pairs

### read amplicon consensus sequence, trim barcode (search for primer sequence and trim flanking sequence)

### pool all reads from each barcode

### assembly



primer_info_file    = primer_pairs






df_barcodes = pd.read_csv(barcode_list, sep='\t', skiprows=0, header=0)
for index, row in df_barcodes.iterrows():
    barcode_name = str(row['f_barcode_name']) + "_" + str(row['r_barcode_name'])



    merged_fasta = open(str(consensus_output) + barcode_name + "/"+ barcode_name + "merged_reads" +".fasta", "w")
    df_primers = pd.read_csv(primer_info_file, sep='\t', skiprows=0, header=0)  ### read primer info
    for index, row in df_primers.iterrows():
        #print row
        f_primer_name, r_primer_name =  str(row['f_primer_name']), str(row['r_primer_name'])
        primer_chr_no, amplicon_start, amplicon_stop = int(f_primer_name.split("_")[1]), int(f_primer_name.split("_")[2]), int(r_primer_name.split("_")[2])
        #seq_path = "./amplicon" + "_" + "_" +  str(amplicon_start) + "_" + str(amplicon_stop) + "/"
        seq_path = str(consensus_output) + barcode_name + "/amplicon_"+ str(primer_chr_no) + "_" +  str(amplicon_start) + "_" + str(amplicon_stop) + "/"
        #print seq_path
        fasta_sequences = SeqIO.parse(open(seq_path + "amplicon_analysis.fasta"), 'fasta') 
        for fasta in fasta_sequences:
            header, sequence = fasta.id.split('_'), str(fasta.seq)
            no_reads = int(header[3][8:])
            if no_reads >= no_reads_threshold:
                merged_fasta.write(">" +str(barcode_name) + "_" +str(primer_chr_no) + "_" +str(amplicon_start)+ "_" +str(amplicon_stop) + "\n")
                merged_fasta.write(sequence[trim_bp:-trim_bp] + "\n")
    merged_fasta.close()
    
    
    
    
    minimus_assembly(barcode_name, amos_path, consensus_output)


'''
### convert fasta to afg format
reads = consensus_output + barcode_name + "/" + barcode_name + "merged_reads"

p = sp.Popen(["%stoAmos" %amos_path,"-s","%s" %reads+".fasta", "-o", "%s" %reads+"_assembly.afg"], stdout=sp.PIPE)
out, err = p.communicate()
### assemble afg reads
p = sp.Popen(["%sminimus" %amos_path, "%s" %reads+"_assembly.afg"], stdout=sp.PIPE)
out, err = p.communicate()
shutil.rmtree(reads+ "_assembly.bnk")
os.remove(reads+ "_assembly.afg.runAmos.log")
os.remove(reads+ "_assembly.afg") 
os.remove(reads+ "_assembly.contig") 
'''

#############################################
#Time to run the code: end timer
############################################################
t1 = time.time()
total = t1-t0
total = ("{0:.2f}".format(round(total,2)))
print "total time to run = ", total, " seconds"


















'''
if __name__ == '__main__':

    #################################################
    ### for each barcode combinations
    required_barcode = 0
    
    
    barcodes_list = [0]
    #################################################
    
    
    

    
    ###Processors to use####
    # if cores == 0 :## Use default 95% of processors as pool
        # nproc = int(multiprocessing.cpu_count()*0.90)
    # else:## As mannually entered by the user
        # nproc = int(cores)
    
    # print("\n#### %s cores reserved for analysis #########" % (str(nproc)))
    # print("#### %s threads assigned to one lib #########\n" % (str(nthread)))
    
    
    
# def worker():
    # """worker function"""
    # print 'Worker'
    # return

# if __name__ == '__main__':
    # jobs = []
    # for i in range(5):
        # p = multiprocessing.Process(target=worker)
        # jobs.append(p)
        # p.start()
    
    
    pooled_locus_read_assembly(loci_info, primer_pairs, barcodes_list)






    
    
    
    ############################################################
    #Time to run the code: end timer
    ############################################################
    t1 = time.time()
    total = t1-t0
    total = ("{0:.2f}".format(round(total,2)))
    print "total time to run = ", total, " seconds"

'''


