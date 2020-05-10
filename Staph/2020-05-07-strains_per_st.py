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

    df = pd.read_csv('/Users/liamcheney/Desktop/Book2.csv')

    infile = open('/Users/liamcheney/Desktop/Untitled.txt').read().splitlines()

    sub = df[df['ST'].isin(infile)]
    sub = sub.sort_values('Contig Count')
    sub.to_csv('/Users/liamcheney/Desktop/UntitledXXX.csv')





if __name__ == '__main__':
    main()
