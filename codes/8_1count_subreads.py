### get # raw reads corresponding to a set of ccs reads
### read alignment accuracy statistics
### Version 0.1.0: 03/19/2017
### Author: Felix Francis (felixfrancier@gmail.com)

############################################################
#Time to run the code: start timer
############################################################
import time
t0 = time.time()

############################################################
#### IMPORT FUNCTIONS
############################################################




############################################################
#### PATHS AND INPUT FILES
############################################################
### raw reads file, paths
path_raw_reads = "/mnt/data27/ffrancis/PacBio_sequence_files/EqAd_DepthCoverage_accuracy_030717/read_nos/"
raw_reads_headers = "m160901_060459_42157_c101086112550000001823264003091775_s1_p0_headers.txt"
raw_reads_file = path_raw_reads + raw_reads_headers

### ccs read id file, paths
ccs_reads_path = "/mnt/data27/ffrancis/PacBio_sequence_files/EqAd_DepthCoverage_accuracy_030717/consensus_from1_50ccs/output_50/"
ccs_reads_filename = "TA_1_25390617_27_F_TA_1_25395472_24_R_reads.txt"
ccs_reads_file = ccs_reads_path + ccs_reads_filename



############################################################
#### FUNCTIONS
############################################################

def count_subreads(raw_reads_file, ccs_id):
	f = open(raw_reads_file)
	contents = f.read()
	f.close()
	# ccs_reads_file
	return contents.count(ccs_id)

def total_subreads_per_ccsidfile(ccs_reads_file, raw_reads_file):
	ccs_file = open(ccs_reads_file)
	raw_read_count = 0
	for ccs_id in ccs_file:
		ccs_id =  ccs_id.strip("\n")
		subreads_count_perccs = count_subreads(raw_reads_file, ccs_id)
		raw_read_count += subreads_count_perccs
	ccs_file.close()
	return raw_read_count


############################################################
#### CODE
############################################################

if __name__ == "__main__":
	print total_subreads_per_ccsidfile(ccs_reads_file, raw_reads_file)



#############################################
#Time to run the code: end timer
############################################################
t1 = time.time()
total = t1-t0
total = ("{0:.2f}".format(round(total,2)))
print "total time to run = ", total, " seconds"


















'''


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


