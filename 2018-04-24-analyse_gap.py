#script aims to iterate over every line, and if each cell has a any locus then select and store this line.
import math
import os
from time import sleep as sl
import pandas as pd
import sys
import glob
import argparse
import matplotlib
matplotlib.use('Agg')
# import matplotlib.pyplot as plt
import shutil
import re

parser = argparse.ArgumentParser(description = "Tool to analyse roary gene_abscence_presence.csv files. Required arguments  [-i][-o][-ra][-s][-c]")
parser.add_argument('-i', '--input', type = str, help = "provide an absolute path to folder containing roary output files.")
parser.add_argument('-ra', '--reference_accession', type = str, help = "provide the accesion number for the reference or genome of interest. Eg GCA_0000001.")
parser.add_argument('-s', '--strain_input', type = int, help = "the number of genomes strains included in roary analysis.")
parser.add_argument('-pa','--previously_analysed', nargs = '?', type = str, help = "provide an absolute path to a file containing a list of genes previously used for in higher conservation schemes.")
parser.add_argument('-o', '--output_files', type = str, help = "provide absoulate path to output")
parser.add_argument('-c', '--convert_cds', type = str, help = "provide a path to CDS-to-gene.txt file. See documentation for further help.")
parser.add_argument('-ssub', '--skip_subfiles', type = str, help = "Use either 'Yes' to skip subfile analysis")

args = parser.parse_args()

strain_input = args.strain_input
reference = args.reference_accession
in_path = args.input
output_files = args.output_files
analysed_ds = args.previously_analysed
cds_to_gene = args.convert_cds
skip_subfiles = args.skip_subfiles

if skip_subfiles == 'Yes':
    pass
else:
    strains_to_keep = math.ceil(int(strain_input) * 0.99) #is the number of strains to include orthologs from. Number is rounded up.
    if os.path.exists(output_files + '/subfiles/'):
        shutil.rmtree(output_files + '/subfiles/')
    os.mkdir(output_files + '/subfiles/')

    if os.path.exists(output_files + '/gene_lists/'):
        shutil.rmtree(output_files + '/gene_lists/')
    os.mkdir(output_files + '/gene_lists/')

    graph_dict = {}
    for filename in sorted(glob.iglob(in_path + '/*'), reverse=True):
        infile = pd.read_csv(filename, index_col=0, header=0, low_memory=False)
        blast_percent = filename[-28:-26]

        sub_df = infile.loc[infile['Avg sequences per isolate'] == 1]
        sub_df = sub_df.loc[sub_df['No. isolates'] >= strains_to_keep]
        sub_df.to_csv(output_files + '/subfiles/subfile_'+str(blast_percent)+'.csv', sep=',')
        graph_dict[blast_percent] = sub_df.shape[0]
        print("Created Sub-dataframe " + filename + " from blast identity " + str(blast_percent) + "%")

    plotting_data = open(output_files + '/plotting_data.tsv', 'w')
    plotting_data.write('Blast Percentage' + '\t' + 'Number of Core Genes' + '\n')
    for a,b in graph_dict.items():
        plotting_data.write(str(a) + '\t' + str(b) + '\n')

    # graph_sort = sorted(graph_dict.items())
    # x,y = zip(*graph_sort)
    # plt.plot(x,y, lw = 2)
    # plt.xlabel('Percentage of Blast Identity')
    # plt.ylabel('Number of Core Genes')
    # plt.title('Impacts of varying blast percentages of core genes numbers')
    # plt.savefig(output_files +'plot.pdf')
    # print('\n' + "Graph was saved as pdf")

#the previously created excel will have all orthologs with no paralogs at the desired cutoff. here the cds locus tags are converted into the gene names provided by a reference input.
already_written = []

for filename in sorted(glob.iglob(output_files + '/subfiles/*'), reverse=True):
    subfile = open(filename, 'r')
    df = pd.read_csv(subfile, index_col=0)
    df = df.fillna("NaN")

    reference_info_input = open(cds_to_gene,'r').readlines()
    ref_dict = {}
    for k in reference_info_input:
        k, v = k.split('\t')
        ref_dict[k] = v.strip('\n')

    blast_percent = filename[-6:-4]
    gene_scheme_out = open(output_files + '/gene_lists/genes_' + blast_percent + '%.txt','w')

    if analysed_ds == None : #will read in a list of previously used genes to prevent duplications
        gene_used_in = []
    else:
        gene_used_in = open(analysed_ds, 'r').read().split('\n') # import genes already added to a previous scheme to prevent duplicating genes across all schemes

    for i in df[reference]:
        if "cds" in i:
            ind = i.index("_") or i.index(":")
            locus = i[:ind]
            if ref_dict[locus] not in gene_used_in and ref_dict[locus] not in already_written:
                gene_scheme_out.write(ref_dict[locus] + '\n')
                already_written.append(ref_dict[locus])
        elif "na" in i:
            pass
    print("Created gene list for " + filename + " from blast identity " + blast_percent + "%")
    gene_scheme_out.close()


gene_used_out = open(output_files + '/used_genes.txt', 'a') #will create a file and add all the genes which were written above. Used for later analysis to prevent duplication genes in schemes.
for i in already_written:
    gene_used_out.write(str(i) + '\n')
gene_used_out.close()
