##Overview of program
##want to create: a MSA for each gene located across the entire input dataset
## 1. iterate through all genomes for each core gene
## 2. extract out the sequence for each core gene from GFF files
## 3. create a MSA for each core gene of that locus from all genomes

##process
## 1. create a dictionary with genome_id[node##:sequence, ...]
## 2. iterate over each core gene in roary output, collect information to find sequence from each dict entry
## 3. save of all sequences for a core gene to a dictionary
## 4. determine which strand orientation the core gene is, and convert (or not) to make same as the reference gene oreintation
## 5. write core gene sequence from all genomes to a single MSA

from time import sleep as sl
import glob
from Bio import SeqIO
import pandas as pd
from Bio.Seq import Seq
import numpy as np

# fnas_input_folder = "/Users/liamcheneyy/Desktop/subset_50_genomes/1_fnas/"
# roary_csv_input = open("/Users/liamcheneyy/Desktop/subset_50_genomes/2_roary_subset.csv", 'r').read().splitlines()
# MSA_out_folder = "/Users/liamcheneyy/Desktop/subset_50_genomes/3_MSA_fna/"
# genomes_subset_list = list([x.strip('.fna') for x in open('/Users/liamcheneyy/Desktop/subset_50_genomes/1_fnas/SNAP_subset_genomes.txt').read().splitlines()])
snap_input_path = "/Users/liamcheneyy/Desktop/8_SNAP_Extra/"

def making_roary_subset(genomes_subset_list, MSA_out_folder):
    df_in = pd.read_csv('/Users/liamcheneyy/Desktop/subset_50_genomes/1_D2C_i96_99_edited_names_percen_para_ana.csv', index_col=False)
    begining_cols = df_in.columns.values.tolist()[0:15]
    df_in = df_in[df_in.columns.intersection(genomes_subset_list + begining_cols)]
    df_in.to_csv(MSA_out_folder + '/roary_subset.csv', index=False)

def fna_data_dict(fnas_input_folder):
    ##creating dict with all fna information
    fna_data = {}
    for filename in glob.iglob(fnas_input_folder + '/*.fna'):
        print(filename)
        genome_id = filename.split('/')[-1][:-4]
        fna_data[genome_id] = {}
        for record in SeqIO.parse(filename, 'fasta'):
            fna_data[genome_id][record.id] = str(record.seq)
    return fna_data

def core_gene_iterator(roary_csv_input, MSA_outfile_path):
    ##iterating over roary csv and creating MSA.fasta file with seq for core gene across all dataset
    genomes_used_df = pd.DataFrame(columns=['Number of Used Genomes', 'Number of Not Used Genomes', 'Number of Excluded Sequences with Ns', 'Number of Excluded Truncated Sequences'])
    fna_data = fna_data_dict(fnas_input_folder)

    for line in roary_csv_input[1:]:
        col = line.split(',')
        genome_used_count = 0
        N_count = 0
        outfilename = col[0]
        print(outfilename)
        MSA_dict = {}
        ref_seq = ''
        for ref_gene in col[15:1949]:
            if 'GCA000006745' in ref_gene:
                ref_genome_id = ref_gene.split('_')[0].strip('"')
                ref_node = ref_gene.split('_')[1]
                ref_start = int(ref_gene.split('_')[3])
                ref_end = int(ref_gene.split('_')[4])
                ref_strand = ref_gene[-2]
                ref_seq = fna_data[ref_genome_id][ref_node][ref_start - 1:ref_end]
                #making all ref seqs in the forward direction
                if ref_strand == '-':
                    ref_seq = Seq(ref_seq).reverse_complement()
                    strand = '+'

                MSA_dict[ref_genome_id + '_' + ref_node] = ref_seq
                # print(ref_genome_id)
                # print(ref_strand, ref_seq)

        for gene in col[15:1949]:
            if '\t' in gene:
                pass

            elif gene != '' and 'GCA000006745' not in gene:
                non_ref_seq = ''
                roary_genome_id = gene.split('_')[0].strip('"')
                start = int(gene.split('_')[3])
                end = int(gene.split('_')[4])
                strand = gene[-2]

                #used to fix bad roary cell strings
                if 'ERR' in roary_genome_id and 'NODE' not in gene.split('_')[1]:
                    roary_node = 'NODE' + gene.split('_')[1]
                else:
                    roary_node = gene.split('_')[1]

                #checking this sequence against the ref sequences
                non_ref_seq = fna_data[roary_genome_id][roary_node][start - 1:end]

                #make all sequnces the same strand
                if strand == '-':
                    non_ref_seq = Seq(non_ref_seq).reverse_complement()

                #removing strains with Ns
                if 'N' in non_ref_seq:
                    N_count = N_count + 1
                else:
                    MSA_dict[roary_genome_id + '_' + roary_node] = str(non_ref_seq)


        len_avg = []
        for key, value in MSA_dict.items():
            len_avg.append(len(value))
        avg_len = (sum(len_avg) / len(len_avg))

        short_count = 0
        remove_keys = []
        for key, value in MSA_dict.items():
            if '000006745' not in key:
                if len(value) >= (avg_len + 3) or len(value) <= (avg_len - 3):
                    remove_keys.append(key)
                    short_count = short_count + 1

        for i in remove_keys:
            del MSA_dict[i]

        #writing MSA outfile
        outfile = open(MSA_outfile_path + outfilename + '.fna', 'w')
        for key, value in MSA_dict.items():
            outfile.write('>' + str(key) + '\n')
            outfile.write(str(value) + '\n')
        outfile.close()

        # write genomes used stats
        genome_used_count = len(MSA_dict.keys())
        genomes_not_used = 1935 - genome_used_count
        genomes_used_df.loc[outfilename] = [genome_used_count, genomes_not_used, N_count, short_count]

    genomes_used_df.to_csv(MSA_outfile_path + 'Genomes_used_stats.csv')

def snap_output_ds(input_path):
    ##not used since changed perl program to increase decimals output. dont need to compute myself
    line_c = 0
    #create list of ds for each file and find average
    for filename in glob.iglob(input_path + '*group_*'):
        core_group = filename.split('/')[-1].strip('summamry.')
        list_of_ds = []
        excluded_lines = 0
        infile = open(filename,'r').read().splitlines()
        for line in infile[1:-2]:
            line_c = line_c + 1
            col = line.split()
            if '-' not in col[-3]:
                ds = float(col[-3])
                list_of_ds.append(ds)
            if '-' in col[-3]:
                excluded_lines = excluded_lines + 1

        print(sum(list_of_ds))
        #calculate the avg ds
        average_ds = np.average(list_of_ds)
        # average_ds = sum(list_of_ds) / len(list_of_ds)
        print(core_group +'\t' + str(round(average_ds,4)))

core_gene_iterator(roary_csv_input, MSA_out_folder)
