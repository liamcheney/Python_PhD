import argparse
from time import sleep as sl
import pandas as pd
import numpy as np
from Bio import SeqIO

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args

def get_first_strain(df,in_genes):
    #select certain rows
    sub_df = df[df.index.isin(in_genes)]
    begining = sub_df.iloc[:,:13]
    genomes = sub_df.iloc[:,13:]
    genomes.dropna(axis=1, inplace=True ) #remove strains not containing gene
    new = pd.merge(begining, genomes, left_index=True, right_index=True)
    new.to_csv('/Users/liamcheneyy/Desktop/sub.csv')

    # # first strain with complete contig
    save_dict = {}
    for strains in list(new.columns.values)[13:]:
        contig_list = []
        strain = new[strains].values.tolist()
        for i in strain:
            contig = i.split('_')[1]
            contig_list.append(contig)

        contig_list = list(set(contig_list))

        save_dict[strains] = len(contig_list)

    print(save_dict)
    strain_1_acc = list(new.columns.values)[13:14][0]
    strain_1 = new[strain_1_acc].values

    save_list = []
    for i in strain_1:
        x = i.split('_')
        save_list.append(x)

    return save_list

def check_contigs(save_list):
    # get list of contigs
    contig_list = [x[1] for x in save_list]
    contig_list = list(set(contig_list))

    if len(contig_list) > 1:
        print('Multiple contigs : ' + str(len(contig_list)))
    else:
        print('Single contig.' + str(len(contig_list)))
def check_sequential(first_strain, genome_in_path):
    # sort contigs based on gene ID
    first_strain.sort(key=lambda x: x[2])

    # first check if all gene IDs are sequential
    genes_ids = [int(x[2]) for x in first_strain]

    true_count = 0
    for element_num in range(0, len(genes_ids) - 1, 1):
        element = genes_ids[element_num]
        next = genes_ids[element_num + 1]
        if (element + 1) == next:
            true_count = true_count + 1

    if true_count == (len(genes_ids) - 1):
        extract_seq(genome_in_path, first_strain)
    else:
        return False
def get_longest_block(first_strain):
    # sort contigs based on gene ID
    first_strain.sort(key=lambda x: x[2])

    #get the longest stretch of DNA
    save = []
    for element_num in range(0, len(first_strain) - 1, 1):

        stretch_save = []
        for row_num in range(0, len(first_strain) - 1, 1):
            row_num = row_num + element_num

            first_id = int(first_strain[row_num][2])
            next_id = int(first_strain[row_num + 1][2])

            # if (first_id + 1) == next_id:
            #     stretch_save.append(first_id)
            #     elif (first_id + 1) != next_id:
            #         save.append(stretch_save)
            #         break





    # for i in save_list:
    #     print(i)

    # #check if genes are next to each other
    # print(strain_1_acc)
    # for row in strain_1:
    #     start = row.split('_')[3]
    #     end = row.split('_')[4]
    #     contig = row.split('_')[1]
    #     gene = row.split('_')[2]
    #     print(contig, gene,start,end,sep='\t')


def extract_seq(genome_in_path, first_strain):
    complete_start = int(first_strain[0][3])
    complete_end = int(first_strain[-1][4])
    contig = first_strain[0][1]
    accession = first_strain[0][0]
    for record in SeqIO.parse(genome_in_path + '/' + accession + '.fna', 'fasta'):
        if contig in record.id:
            record.seq = (record.seq)[complete_start:complete_end]
            SeqIO.write(record,genome_in_path + '/' + contig + '_' + str(complete_start) + '_' + str(complete_end) + '.fasta','fasta')
def main():
    args = parseargs()

    #read in roary csv with extra info
    df = pd.read_csv('/Users/liamcheneyy/Desktop/gap_with_details.csv', index_col=0, low_memory=False)

    #read in list of interest genes
    in_genes = open('/Users/liamcheneyy/Desktop/genes.txt').read().splitlines()
    in_genes = [x for x in in_genes]

    genome_in_path = '/Users/liamcheneyy/Desktop/all_fasta/'

    #find first strain with gene on a single contig
    first_strain = get_first_strain(df, in_genes)

    #check contigs, if cant find strain with single contig
    check_contigs(first_strain)

    # check if sequential
    sequential_check = check_sequential(first_strain,genome_in_path)

    #if not sequential then get longest run of genes possible
    # get_longest_block(first_strain)

    #TODO get the strains which has the most amount of same contigs for genes. get values.
    #TODO return this strain for further analysis

    #TODO get longest stretch of alignment and return

if __name__ == '__main__':
    main()

#logic: finding blocks of genes
#take a list of genes
#create matrix with extra information added
#use positions from first genome, check if are next to each other.


##for getting sequence
    #check if genome exists, where they are all co-linear
        #check if they on same contig and next to each other
