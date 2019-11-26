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
    df = pd.read_csv('/Users/liamcheneyy/Desktop/vcseventh_22/grapetree/seventh/MGT_isolate_data.txt', sep='\t', low_memory=False, index_col=0)

    out_path = '/Users/liamcheneyy/Desktop/new_combinations.tsv'

    #create combinations of alleles
    ctxB = ['ctxB1','ctxB3','ctxB4','ctxB5','ctxB7']
    tcpA = ['tcpA_A226',	'tcpA_el', 'tcpA_cla']
    rstR = ['rstR_el', 'rstR_cla','rst_cla/el']
    combination = list(itertools.product(ctxB, tcpA, rstR))
    combinations = [list(x) for x in combination]

    #prelim check. get numbers of each allele combination
    freq_dict = {}
    comb_count = 1
    keep_combination = []
    for comb in combinations:
        sub_df = df[comb]

        for i in comb:
            sub_df = sub_df[sub_df[i] == True]

        strains_num = sub_df.shape[0]

        if strains_num > 0:
            comb_string = '-'.join(comb)
            freq_dict[comb_string] = {'comb':comb, 'size':strains_num, 'strains':list(sub_df.index)}
            keep_combination.append(comb_string)

        comb_count = comb_count +1

    zip(keep_combination)

    #get shorthand list results
    results = []
    for key, value in freq_dict.items():
        string = '-'.join(value['comb'])
        num = value['size']
        results.append(num)


    #create df and fill TRUE or FALSE for each strain combiniation
    save_dict = {}
    len_of_combinations = len(keep_combination)
    for strain in list(df.index):
        save_dict[strain] = []
        for comb in freq_dict.keys():
            if strain in freq_dict[comb]['strains']:
                save_dict[strain].append('TRUE')
            elif strain not in freq_dict[comb]['strains']:
                save_dict[strain].append('FALSE')


    with open(out_path,'w') as out:
        out.write('ID' + '\t')
        for i in keep_combination:
            out.write(i + '\t')
        out.write('\n')

        for key, value in save_dict.items():
            out.write(key + '\t')
            for k in value:
                out.write(k + '\t')
            out.write('\n')




if __name__ == '__main__':
    main()
