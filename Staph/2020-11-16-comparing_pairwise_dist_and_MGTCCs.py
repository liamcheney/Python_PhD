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

def find_same_CC_and_ODC(pairwise_dict, CC_dict):

    ##find same CC and PWs
    print(CC_dict['62'])
    #find majority matches
    for CC, CC_list in CC_dict.items():
        CC_strains_num = len(CC_list)
        save_list = []
        # if CC not in complete_list:
        if CC == '62':
            for ODC, ODC_list in pairwise_dict.items():

                #find strains in both
                for element in ODC_list:
                    if element in CC_list:
                        print(element)
                        print(ODC_list, CC_list)

                        save_list.append(element)

            # percen1 = str(len(save_list)/CC_strains_num * 100)
            # print(CC, percen1, CC_list, ODC_list)



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
