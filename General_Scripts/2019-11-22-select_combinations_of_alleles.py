import argparse
from time import sleep as sl
import itertools
import pandas as pd
import sys

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args


def main():
    args = parseargs()

    #read in dataframe
    df = pd.read_csv('/Users/liamcheneyy/Desktop/MGT_isolate_data.txt', sep='\t', low_memory=False, index_col=0)

    #create combinations of alleles
    ctxB = ['ctxB1','ctxB3','ctxB4','ctxB5','ctxB7']
    tcpA = ['tcpA_el_A226',	'tcpA_el_WT', 'tcpA_cla_WT']
    rstR = ['rstR_el', 'rstR_cla']
    combination = list(itertools.product(ctxB, tcpA, rstR))
    combinations = [list(x) for x in combination]

    #prelim check. get numbers of each allele combination
    freq_dict = {}
    for comb in combinations:
        sub_df = df[comb]

        for i in comb:
            sub_df = sub_df[sub_df[i] == True]
        print(comb,sub_df.shape)
        sl(1)


        # print(true_sub_df.shape)
        # true_sub_df = true_sub_df[true_sub_df[comb[0] == 'TRUE']]
        # print(true_sub_df)
        # strains_num = true_sub_df.shape
        # print(comb, strains_num)
        # sl(1)

    print(combinations)

if __name__ == '__main__':
    main()
