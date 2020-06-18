import argparse
from time import sleep as sl
from Bio import SeqIO
import glob
import pandas as pd
import random
from pycountry_convert import country_alpha2_to_continent_code, country_name_to_country_alpha2

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args


def main():
    args = parseargs()


    infile = open('/Users/liamcheneyy/Desktop/meta_7-gene-alleles.txt').read().splitlines()

    for line in infile[1:]:
        col = line.split('\t')
        code = col[6].strip('"')
        if code.isalpha():
            print(country_alpha2_to_continent_code(code))
        else:
            print()



if __name__ == '__main__':
    main()