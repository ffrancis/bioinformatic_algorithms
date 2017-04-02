### randomly sample ccs ids from a file of all ccs ids for a given amplicon
### only for 21-25 th replicates..so remove make directory options and only add these new replicates to the same original directories
### this will be done for a range of # of ccs ids. example 1-35
### Version 0.1.0: 032417
### Author: Felix Francis (felixfrancier@gmail.com)

############################################################
#Time to run the code: start timer
############################################################
import time
t0 = time.time()

############################################################
#### IMPORT FUNCTIONS
############################################################

# from Bio import SeqIO
import pandas as pd
import numpy as np
import os


############################################################
#### PATHS AND INPUT FILES
############################################################

### path to all ccs id files
all_ccs_path = "/mnt/data27/ffrancis/PacBio_sequence_files/EqAd_DepthCoverage_accuracy_030717/consensus_from1_50ccs_with_20sampliing/output_all/"


primer_pairs_file    = "primer_pairs_info.txt"

fofn_pacbio_raw_reads   = "/mnt/data27/ffrancis/PacBio_sequence_files/EqPCR_raw/F03_1/Analysis_Results/m160901_060459_42157_c101086112550000001823264003091775_s1_p0.bas.h5"
processors  = 12
walltime    = 400
no_replicates = 20
maxno_ccs_ids = 40


############################################################
#### FUNCTIONS
############################################################

# function for sampling with replacement
def random_no_sampling(no_ccs_reads, no_samples_required):
	return np.random.choice(no_ccs_reads, no_samples_required)



############################################################
#### CODE
############################################################


# os.makedirs("./output_sampling/")
os.makedirs("./output_sampling_1_20_38/")


for i in xrange(1, maxno_ccs_ids+1):
	no_ccs_ids = i
	
	# no_ccs_ids = 2			### modify later to iterate over 1 - n number of no_ccs_ids

	output_directory = "./output_sampling_1_20_38/output_" + str(no_ccs_ids) +"/"
	os.makedirs(output_directory)

	### read primer info file
	df_primers = pd.read_csv(primer_pairs_file, sep='\t', skiprows=0, header=0)  ### read primer info
	# for index, row in df_primers.iterrows():
		# f_primer_name, r_primer_name =  str(row['f_primer_name']), str(row['r_primer_name'])
		# amplicon_name = str(f_primer_name) + "_" + str(r_primer_name)
		### filename = f_primer_name + "_" + r_primer_name + "_reads.txt"
		# filename = amplicon_name + "_reads.txt"





	with open("./output_sampling_1_20_38/" +  "consensus_calling_"+ str(no_ccs_ids) + ".sh", "w") as shell_script_output: 
		cwd = os.getcwd()
		shell_script_output.write("#!/bin/sh" + "\n")

		shell_script_output.write("#PBS -N consensus_calling_"+ str(no_ccs_ids) +"\n")
		shell_script_output.write("#PBS -r n"+ "\n")
		shell_script_output.write("#PBS -l walltime="+ str(walltime) +":00:00" + "\n")

		shell_script_output.write("#PBS -l nodes=1:ppn=" + str(processors) + "\n")

		shell_script_output.write("#PBS -d " + cwd + "./output_sampling_1_20_38/" + "\n"+  "\n")
		shell_script_output.write("export SMRT=/opt/smrtanalysis" +  "\n")
		shell_script_output.write("export INPUT=" + fofn_pacbio_raw_reads +  "\n")
		shell_script_output.write("alias smrtwrap='$SMRT/smrtcmds/bin/smrtwrap'" + "\n" + "\n")
		# output_dir = "./output_" + str(no_ccs_ids) +"/"

		
		### iterate per amplicon
		for index, row in df_primers.iterrows():
			f_primer_name, r_primer_name =  str(row['f_primer_name']), str(row['r_primer_name'])
			amplicon_name = str(f_primer_name) + "_" + str(r_primer_name)
			chr_no, amplicon_start, amplicon_stop =  f_primer_name.split("_")[1], f_primer_name.split("_")[2], r_primer_name.split("_")[2]
			filename = amplicon_name + "_reads.txt"

			### read file with all ccs reads
			f=open(all_ccs_path + filename)
			lines=f.readlines()
			no_ccs_reads = len(lines)
			if len(lines) > 0:
			

				### iterate per replicate
				for rep in xrange(1, no_replicates+1):
					reads_filename = amplicon_name + "_" + str(rep) + "_reads.txt"
					
					### command for consensus calling
					shell_script_output.write("### laa error corrected consensus calling "+ amplicon_name + "_ reads" +  "\n")
					shell_script_output.write("smrtwrap ConsensusTools.sh AmpliconAnalysis $INPUT --noPhasing --noClustering -n " + str(processors) + " --whiteList " + "./output_" + str(no_ccs_ids) +"/" + reads_filename + " -o " + "./output_" + str(no_ccs_ids) +"/" + "amplicon_" + str(chr_no) + "_" + str(amplicon_start) + "_" + str(amplicon_stop) + "_" + str(rep) + "/" + "\n" + "\n")
					
					### print randomly sampled (with replacement) ccs read ids
					with open(output_directory + reads_filename, "w") as output:
						list_line_nos =  random_no_sampling(no_ccs_reads, no_ccs_ids)
						for i in list_line_nos:
							output.write(lines[i])
			


	
	
	
	
	
	
	
#############################################
#Time to run the code: end timer
############################################################
t1 = time.time()
total = t1-t0
total = ("{0:.2f}".format(round(total,2)))
print "total time to run = ", total, " seconds"




