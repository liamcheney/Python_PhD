import argparse
from time import sleep as sl
import pandas as pd

def SNP_to_gene_compare(gene_ref_pos_df, snp_df):

    for position in snp_df['Position']:

        #check what gene positions is in
        sb = ((gene_ref_pos_df[gene_ref_pos_df['Start'] < position]) & (gene_ref_pos_df[gene_ref_pos_df['End'] > position]))
        print(sb)


def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args


def main():
    args = parseargs()

    #input of genes
    genes_positions_path = '/Users/liamcheneyy/Desktop/ref_sample.txt'
    gene_ref_pos_df = pd.read_csv(genes_positions_path, sep='\t')

    #branch number for da long script
    branch_num = 7

    # input SNPS tsv
    snp_input_path = '/Users/liamcheneyy/Desktop/snp_samples.txt'
    snp_df = pd.read_csv(snp_input_path, sep='\t')
    snp_df = snp_df[snp_df['Branch'] == branch_num]

    #compare the da long SNPs and my script alleles
    SNP_to_gene_compare(gene_ref_pos_df, snp_df)

    #compare genes to ones from input



if __name__ == '__main__':
    main()
