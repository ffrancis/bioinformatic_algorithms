#compute_gene_space

#path        ="./"
#input_file  ="sample.txt"


path        ="/net/biohen27/data/ffrancis/ZmDB_RefGenv3_annotation/"
input_file  ="zea_mays.protein_coding.gff"

output      = open("gene_space_info.txt", "w")


chr_dict    =   {}
total_gene_space    =   0

with open (path + input_file, "r") as input_file:
    lines       =   input_file.readlines()
    for line in lines:
        if not (line.startswith('#') or line.startswith('scaffold')):
            cols    =   line.split()
            chr_no, annotation, start, stop  =   cols[0], cols[2], cols[3], cols[4]
            if annotation   ==   "gene":
                gene_len    = (int(stop)-int(start))+1
                total_gene_space    += gene_len
                if chr_no in chr_dict:
                    chr_dict[chr_no] = (chr_dict[chr_no] + gene_len)
                else:
                    chr_dict[chr_no] = gene_len

#print chr_dict

output.write ("Total gene space    = " + str(total_gene_space) + " bp" + "\n" + "\n")


output.write ("chr_no" + "\t" + "gene_space(bp)"+ "\n")


for chr_no, gene_space in chr_dict.iteritems():
    #print chr_no, gene_space
    output.write (str(chr_no) +"\t" + str(gene_space) + "\n")
                   
                   
                   
output.close

