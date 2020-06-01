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

    df = pd.read_csv('/Users/liamcheney/Desktop/MGT9_allele_profiles.txt',sep='\t')

    for col in df:
        if "#" not in col:
            sub = df[col].astype(int)
            uniques = sub.value_counts().shape[0]

            # for i,r in uniques.iteritems():
            #     print(i,r)
            print(col, uniques, sep='\t')

if __name__ == '__main__':
    main()