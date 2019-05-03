import argparse
from time import sleep as sl
import subprocess
from Bio.Blast.Applications import NcbiblastnCommandline


def cerate_files_dicts():


def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-db", "--database_name", required=True,
                        help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args

def main():

    args = parseargs()

if __name__ == '__main__':
    main()

#take in the files for both ctxb alleles and the query genome
#build a blast database

