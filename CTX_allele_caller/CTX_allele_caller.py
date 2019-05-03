import argparse
from os import path
import sys
from time import sleep as sl
import subprocess
from Bio.Blast.Applications import NcbiblastnCommandline


def parseargs(set_wd):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-i", "--input_genome",
                        help="Path to a single input genome (nucleotide) for analysis.")
    parser.add_argument("-o", "--output_folder", default=set_wd + "/output/",
                        help="sql database to search (eg. vcseventh)")
    parser.add_argument("-f", "--input_genomes_file",
                        help="A comma separated list of genomes (nucleotide) to analyse.")

    args = parser.parse_args()

    return args

def blast_input_against_db():

    print()

def main():

    set_wd = path.dirname(path.abspath(__file__))
    args = parseargs(set_wd)



if __name__ == '__main__':
    main()


#blast against database
#get results

