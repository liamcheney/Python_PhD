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
import  glob
from Bio import SeqIO
import pandas as pd
from Bio.Seq import Seq

fnas_input_folder = "/Volumes/Liam's HDD/PhD/2_HGT_seventh_pandemic/2019-01-30-core_genes_DNDS/fnas"
roary_csv_input = open("/Users/liamcheneyy/Desktop/D2C_i96_99_edited_names_percen_para_ana.csv", 'r').read().splitlines()
MSA_out_folder = "/Users/liamcheneyy/Desktop/untitled/"


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
        # print(outfilename)
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

                ## aligner reference sequence against other sequences
                # align_count = 0
                # for x, y in zip(ref_seq[:100], non_ref_seq[:100]):
                #     if x == y:
                #         align_count = align_count + 1
                #
                # if align_count <= (0.70 * 100):
                #     non_ref_seq = Seq(non_ref_seq).reverse_complement()

        #remove truncated or shorter sequence
        # len_list = []
        # for key, value in MSA_dict.items():
        #     len_list.append(len(value))
        # avg_len = round(sum(len_list) / len(len_list),0)
        #
        # short_count = 0
        # remove_keys = []
        # for key, value in MSA_dict.items():
        #     if len(value) < (0. * avg_len):
        #         remove_keys.append(key)
        #         short_count =  short_count + 1
        #
        # for i in remove_keys:
        #     del MSA_dict[i]

        len_avg = []
        for key, value in MSA_dict.items():
            len_avg.append(len(value))
        avg_len = (sum(len_avg) / len(len_avg))

        short_count = 0
        remove_keys = []
        for key, value in MSA_dict.items():
            if 'GCA_00006745'
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

core_gene_iterator(roary_csv_input, MSA_out_folder)

##below is older method which created temp files and then  read in and continued. final product had genes of different lengths
## couple of bugs not figured out.
# genomes_folder = '/Users/liamcheneyy/Desktop/2018-04-02-Gffs_D2C/fnas/'
# roary_input_path = '/Users/liamcheneyy/Google Drive/Honours/4All_Core_Genome/D2C_i96_99_no_ref_blanks_percen_para_ana.csv'
# ref_table_path = '/Users/liamcheneyy/Google Drive/Honours/Reference Info/GCF_000006745.gff'
# ref_genes_path= '/Users/liamcheneyy/Desktop/genes.txt'
# fixing_folder_path = '/Users/liamcheneyy/Desktop/2018-04-02-Gffs_D2C/fixing/'
#
# MSA_out_folder = '/Users/liamcheneyy/Desktop/2018-04-02-Gffs_D2C/MSA_out/'
# fixed_MSA_out_path ='/Users/liamcheneyy/Desktop/2018-04-02-Gffs_D2C/fixed_MSA_out/'
# msa_input_path = '/Users/liamcheneyy/Desktop/2018-04-02-Gffs_D2C/MSA_out/'
#
#
# def split_fna_and_gff(genomes_folder):
#     for filename in glob.iglob(genomes_folder + '/*.gff'):
#         file = open(filename, 'r').read()
#         parts = file.split('##FASTA' + '\n')
#
#         outpath = filename[0:-4] + '.fna'
#         outfile = open(outpath, 'w')
#         outfile.write(parts[1])
#         print(filename)
#         outfile.close()
# # split_fna_and_gff(genomes_folder)
#
# def fix_roary_ref_core_genes(ref_table, genes_list_path):
#     ref_dict = {}
#     infile = open(ref_table, 'r').read().splitlines()
#     for line in infile:
#         if "RefSeq" in line:
#             col = line.split('\t')
#             if 'CDS' in col[2]:
#                 cds = col[8].split(';')[0].split('__')[0][6:]
#                 chromosome = col[0]
#                 start = col[3]
#                 end = col[4]
#                 accession = 'GCA_000006745'
#                 new = accession + '_' + chromosome + '_' + cds + '_' + start + '_' + end
#                 ref_dict[cds] = new
#
#     ref_in_list = list(open(genes_list_path, 'r').read().splitlines())
#     for line in ref_in_list:
#         col = line.split('_')[0][3:]
#         if col in ref_dict.keys():
#             print(ref_dict[col])
#         else:
#             print()
# # fix_roary_ref_core_genes(ref_table_path, ref_genes_path)
#
# def fix_length_fnas(fixing_folder_path):
#     for filename in glob.iglob(fixing_folder_path + '/*fna'):
#         seq_dict = {}
#         accession = filename.split('/')[-1].rstrip('.fna')
#         for record in SeqIO.parse(filename, 'fasta'):
#             node = record.id.split('_')[0] + '_' + record.id.split('_')[1]
#             record.id = accession + '_' + node
#             record.name = ''
#             seq_dict[record.id] = record.seq
#
#         outfile = open(filename[:-20] + 'fixed/' + accession + '.fna','w')
#         for key, value in seq_dict.items():
#             outfile.write('>' + key + '\n')
#             outfile.write(str(value) + '\n')
# # fix_length_fnas(fixing_folder_path)
#
# def gene_finder(genomes_path, roary_input_path, MSA_out_path):
#
#     #create dict with the nodes (scaffolds) for each genome
#     genomes_dict = {}
#     print('Building genomes dictionary.')
#
#     for genome in glob.iglob(genomes_path + '*.fna'):
#         genome_accession = genome.split('/')[-1].rstrip('.fna')
#         genomes_dict[genome_accession] = {}
#         for record in SeqIO.parse(genome, 'fasta'):
#             genomes_dict[genome_accession][record.id] = record.seq
#
#     #iterate over roary output and find sequences in dict
#     roary_in = open(roary_input_path, 'r').read().splitlines()
#     line_count = 1
#     for line in roary_in[1:]:
#         print('Line count = ' + str(line_count))
#         line_count = line_count + 1
#         filename_out = 'none' # created MSA are named after ref gene
#         genes_seq_dict = {} # dict to save all sequences for each gene
#
#         excluded_genes_dict = {}
#         excluded_genes = 0
#         excluded_genes_list = []
#
#         col = line.split(',')
#         total_genomes = len(col[15:-9])
#         genome_count = 0
#         for i in col[15:-9]:
#             genome_count = genome_count + 1
#
#             if 'ERR' in i and '\t' not in i:
#                 accession = i.split('_')[0]
#                 node = accession + '_' + i.split('_')[1] + '_' + i.split('_')[2]
#                 locus_tag = accession + "_" + i.split('_')[3]
#                 start_bp = int(i.split('_')[4])
#                 end_bp = int(i.split('_')[5])
#
#                 if 'length' in i:
#                     node = accession + '_' + 'NODE' + '_' + i.split('_')[1]
#                     locus_tag = accession + "_" + i.split('_')[3]
#                     start_bp = int(i.split('_')[4])
#                     end_bp = int(i.split('_')[5])
#
#                 core_gene_seq = genomes_dict[accession][node][start_bp:end_bp + 1]
#                 genes_seq_dict[locus_tag] = core_gene_seq
#
#             if 'GCA' in i and '\t' not in i:
#                 accession = i.split('_')[0] + '_' + i.split('_')[1]
#                 node = i.split('_')[2] + '_' + i.split('_')[3]
#                 locus_tag = accession + "_" + i.split('_')[4]
#                 start_bp = int(i.split('_')[5])
#                 end_bp = int(i.split('_')[6])
#
#                 if 'GCA_000006745' in i:
#                     filename_out = locus_tag
#                     print(filename_out)
#
#                 if len(i.split('_')) == 7:
#                     accession = i.split('_')[0] + '_' + i.split('_')[1]
#                     node = i.split('_')[2]
#                     locus_tag = accession + "_" + i.split('_')[3]
#                     start_bp = int(i.split('_')[4])
#                     end_bp = int(i.split('_')[5])
#
#                 core_gene_seq = genomes_dict[accession][node][start_bp:end_bp + 1]
#                 genes_seq_dict[locus_tag] = core_gene_seq
#
#
#             if '\t' in i:
#                 excluded_genes = excluded_genes + 1
#                 excluded_genes_list.append(str(i))
#
#             # percen_genome_done = round((genome_count / total_genomes) * 100, 0)
#             # print('Genomes finshed  = ' + str(percen_genome_done) + '%.')
#
#         excluded_genes_dict['number of genes excluded'] = excluded_genes
#         excluded_genes_dict['genes'] = excluded_genes_list
#
#         excluded_genes_out = open(MSA_out_path + filename_out + '_not_used_genes.txt', 'w')
#         excluded_genes_out.write('Not used genes = ' + str(excluded_genes_dict['number of genes excluded']) + '\n')
#         for i in excluded_genes_dict['genes']:
#             excluded_genes_out. write((i + '\n'))
#         excluded_genes_out.close()
#
#         MSA_out = open(MSA_out_path + filename_out + '.fna', 'w')
#         print(MSA_out_path + filename_out + '.fna')
#         for key, value in genes_seq_dict.items():
#             MSA_out.write('>' + key + '\n')
#             MSA_out.write(str(value) + '\n')
#         MSA_out.close()
# # gene_finder(genomes_folder, roary_input_path, MSA_out_folder)
#
# def clean_genes_MSA(MSA_input, fixed_MSA_output_path):
#     start_codon = 'ATG'
#     stop_codons = ['TAA', 'TGA', 'TAG']
#
#     for filename in glob.iglob(MSA_input + '/*.fna'):
#         accession_locus = filename.split('/')[-1][:-4]
#         outfile = open(fixed_MSA_output_path + accession_locus + '_' + 'fixed.fna', 'w')
#         seq_dict = {}
#         for record in SeqIO.parse(filename, 'fasta'):
#             forward_strand = True
#             #fixing start position
#             for i in range(0,3,1):
#                 codon_f = str(record.seq[i:i+3])
#                 codon_r = str(record.seq.reverse_complement()[i:i+3])
#                 if codon_f in start_codon:
#                     new_seq = record.seq[i:]
#                     record.seq = new_seq
#                     forward_strand = True
#
#                 if codon_r in start_codon:
#                     new_seq = record.seq.reverse_complement()[i:]
#                     record.seq = new_seq.reverse_complement()
#                     forward_strand = False
#
#             #fixing end position
#             final_codon_inx = 0
#             for i in range(0, len(record.seq), 1):
#                 if forward_strand == True:
#                     codon_f = str(record.seq[i:i + 3])
#                     if codon_f in stop_codons:
#                         final_codon_inx = i
#
#                 if forward_strand == False:
#                     codon_r = str(record.seq.reverse_complement()[i:i + 3])
#                     if codon_r in stop_codons:
#                         final_codon_inx = i
#
#             fixed_seq = str(record.seq[:final_codon_inx+3])
#             seq_dict[record.id] = fixed_seq
#
#         #writing all sequences to outfile
#         for key, value in seq_dict.items():
#             outfile.write('>' + str(key) + '\n')
#             outfile.write(value + '\n')
#         outfile.close()
#
# split_fna_and_gff(genomes_folder)
# # clean_genes_MSA(msa_input_path, fixed_MSA_out_path)
