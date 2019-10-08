import argparse
from time import sleep as sl
from Bio import  SeqIO
import pandas as pd

def  create_alleles(all_alleles_path, wanted_genes, out_path):

    want_genes = open(wanted_genes,'r').read().splitlines()

    save_list = []
    input = SeqIO.parse(all_alleles_path, 'fasta')
    for record in input:
        name = record.id.split(':')[0]
        if name in want_genes:
            save_list.append(record)

    SeqIO.write(save_list, out_path + "/mgt234_ref_alleles.fna", "fasta")

def create_loci_locations(all_alleles_path, wanted_genes, out_path):
    want_genes = open(wanted_genes,'r').read().splitlines()

    df = pd.read_csv(all_alleles_path, sep='\t')
    # sub_df = df.loc[want_genes]
    # sub_df.to_csv(out_path + '/lociLocationsInRef.csv')
    sub_df = df[df['Locus'].isin(want_genes)]
    sub_df.to_csv(out_path + '/mgt234_lociLocationsInRef.txt', sep='\t', header=False, index=False)

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args


def main():
    args = parseargs()

    ##reading in files
    all_alleles_path = '/Users/liamcheneyy/Desktop/allele_testing/species_ref_alleles.fna'
    all_locations_path = '/Users/liamcheneyy/Desktop/allele_testing/species_lociLocationsInRef.txt'
    wanted_genes = '/Users/liamcheneyy/Desktop/allele_testing/wanted.txt'
    out_path = '/Users/liamcheneyy/Desktop/allele_testing/'

    #create alleles fasta
    create_alleles(all_alleles_path, wanted_genes, out_path)

    #create locations
    create_loci_locations(all_locations_path, wanted_genes, out_path)

if __name__ == '__main__':
    main()
