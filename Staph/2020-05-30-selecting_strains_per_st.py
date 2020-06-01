import argparse
from time import sleep as sl
from Bio import SeqIO
import glob
import pandas as pd
import random

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args


def main():
    args = parseargs()

    STs_in_df = open('/Users/liamcheney/Desktop/strains_st_nm.txt').read().splitlines()

    df = pd.read_csv('/Users/liamcheney/Desktop/meta_7-gene-alleles.txt',sep='\t')

    save_list = []
    for line in STs_in_df:
        col = line.split('\t')
        ST = col[0]
        size = int(col[1])

        sub_df = df[df['ST'] == ST]
        sample = sub_df.sample(n = size)

        save_list.append([ST, sample['ID'].tolist()])

    with open('/Users/liamcheney/Desktop/genomes_STs.txt','w') as out:
        for element in save_list:
            print(element[0])
            for item in element[1]:
                out.write(element[0] + '\t' + item + '\n')


if __name__ == '__main__':
    main()