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

def make_dict(CC_in):
    CC_dict_in = {}
    for line in CC_in[1:]:
        col = line.split('\t')
        strain = col[0]
        CC = col[-1]
        if CC not in CC_dict_in.keys():
            CC_dict_in[CC] = [strain]
        else:
            CC_dict_in[CC].append(strain)

    return CC_dict_in

def find_exact_same(pairwise_dict, CC_dict):
    print(f'Finding exact matches.')
    exact_CC_list = []
    exact_ODC_list = []

    for CC, CC_list in CC_dict.items():
        for ODC, ODC_list in pairwise_dict.items():
            if set(CC_list) == set(ODC_list):
                # print(CC, ODC, len(ODC_list), sep='\t')
                exact_CC_list.append(CC)
                exact_ODC_list.append(ODC)

    for rm_el in exact_CC_list:
        del CC_dict[rm_el]

    for rm_ODC in exact_ODC_list:
        del pairwise_dict[rm_ODC]

    return pairwise_dict, CC_dict

def find_threshold_same(pairwise_dict, CC_dict):
    print(f'Finding threshold matches.')
    threshold = 0.90
    #idea, find the best percentage to select by calculating them all, and then choose best distribution
    ##check

    for CC, CC_list in CC_dict.items():
        overlap_list = []
        total_strains=len(CC_list)
        if CC == '3493':
            for ODC, ODC_list in pairwise_dict.items():
                if len(ODC_list) <= len(CC_list):
                    for strain in ODC_list:
                        if strain in CC_list:
                            overlap_list.append(strain)

                overlap_percent = len(overlap_list) / total_strains * 100
                print(CC, overlap_percent, len(CC_list), len(ODC_list), overlap_list, CC_list, ODC_list)

    return pairwise_dict, CC_dict

def find_same_CC_and_ODC(pairwise_dict, CC_dict):

    ##find exact same CC and ODC clusters
    pairwise_dict, CC_dict = find_exact_same(pairwise_dict, CC_dict)

    ##find similar with X% overlap
    pairwise_dict, CC_dict = find_threshold_same(pairwise_dict, CC_dict)

def main():
    args = parseargs()


    ##Pairwise distances inputs
    PW_in = open('/Users/liamcheney/Desktop/MGT5_CC_clusters_saureus_4_ODC1_2.txt').read().splitlines()
    pairwise_dict = make_dict(PW_in)

    ##MGT CC in
    CC_in = open('/Users/liamcheney/Desktop/MGT_isolate_data.txt').read().splitlines()
    CC_dict = make_dict(CC_in)

    find_same_CC_and_ODC(pairwise_dict, CC_dict)

if __name__ == '__main__':
    main()
