import argparse
from Bio import SeqIO
from time import sleep as sl
from Bio.Blast.Applications import NcbiblastnCommandline
import glob


def iterator(query_folder, subject_path):
    save_dict = {}
    for filename in glob.iglob(query_folder + "/*"):
        locus = filename.split('/')[-1].split('.')[0]
        blast_result = blaster(filename,subject_path)
        keep_result = blast_result[0]
        save_dict[locus] = keep_result

        print(keep_result)

def blaster(filename,subject_path):

    # blast query database against genome
    blast_string = NcbiblastnCommandline(task="blastn", query=filename, subject=subject_path, outfmt=6)
    out, err = blast_string()
    out = out.splitlines()

    return out


def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args
def main():
    args = parseargs()

    #made for comparing two multi fastas against each other. will iterate over query and blast against subject and return top result.

    #input_against_database
    query_folder = "/Users/liamcheneyy/Desktop/cgMLST_loci/"
    subject_path = "/Users/liamcheneyy/Desktop/my_alleles_cat.fasta"

    #read over input sequences and blast against another input file
    iterator(query_folder, subject_path)


if __name__ == '__main__':
    main()
