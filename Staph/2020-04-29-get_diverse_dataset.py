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

    df = pd.read_csv('/Users/liamcheneyy/Desktop/meta_7-gene-alleles.txt', sep='\t')

    ##will create smaller DFs of certain STs, then sorted based on metadata, contigs, year and country, then select.

    ##get list of unique STs
    ST_list = df['ST'].unique().tolist()

    save_list = []
    ##iterate over each ST
    for ST in ST_list:
        sub_df = df[df['ST'] == ST]
        sub_df = sub_df.sort_values(by=['Contig Number', 'Year', 'Country'])
        dropped = sub_df.dropna()

        #if has both year and country metadata, choose lowest contig
        if dropped.shape[0] > 0:
            keep = dropped.head(1)
            accession = keep.values.tolist()
            save_list.append(accession)

        #partial metadata, choose lowest contig with some metadata

        #no metadata, choose lowest contig
        else:
            keep = sub_df.head(1)
            accession = keep.values.tolist()
            save_list.append(accession)

    for i in save_list:
        for q in i:
            print(*q)

    # df = pd.read_csv('/Users/liamcheneyy/Desktop/meta_7-gene-alleles.txt', sep='\t')
    #
    # infile = list(open('/Users/liamcheneyy/Desktop/Untitled.txt').read().splitlines())
    #
    # sub = df[df['ST'].isin(infile)]
    # sub.to_csv('/Users/liamcheneyy/Desktop/XX.csv')
    # print(infile)
    #



if __name__ == '__main__':
    main()
