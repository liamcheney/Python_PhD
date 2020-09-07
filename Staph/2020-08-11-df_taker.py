import argparse
from time import sleep as sl
import pandas as pd
import glob

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args


def main():
    args = parseargs()

    input_path = "/Users/liamcheneyy/Desktop/Preferences.xlsx"
    info = pd.read_excel(open(input_path, 'rb'), index_col=0, sheet_name='Preferences').fillna("none")


    for filename in glob.iglob('/Users/liamcheneyy/Desktop/loci/MGT*_gene_accessions.txt'):
        level = filename.split('/')[-1].split('_')[0]
        MGT_df = pd.read_csv(filename, sep='\t', header=None)
        wants = ['Allele_Changes','Length']
        loci = list(MGT_df.iloc[:,0])
        sub_df = info[info.index.isin(loci)]
        sub_df = sub_df[wants]
        sub_df.index.names = [level]
        sub_df.to_csv('/Users/liamcheneyy/Desktop/loci/both_attr_' + str(level) + '.txt', sep='\t')
        print(level)



if __name__ == '__main__':
    main()
