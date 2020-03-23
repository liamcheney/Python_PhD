import argparse
from time import sleep as sl
import pandas as pd

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args


def main():
    args = parseargs()

    infile_path = "/Users/liamcheneyy/Desktop/vcseventh_26/grapetree/seventh/MGT8_allele_profiles.tsv"
    df = pd.read_csv(infile_path,sep='\t')

    for column in df:
        if "#" not in column:
            sub = df[column]
            count = sub.nunique()
            print(column, count,sep='\t')

if __name__ == '__main__':
    main()
