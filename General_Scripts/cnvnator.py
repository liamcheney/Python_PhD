import argparse
from time import sleep as sl
from subprocess import check_call
import sys

def calc_cnv_bin_size(args):

    strain = ''
    if '/' in args.input:
        strain = args.input.split('/')[-1].split('.')[0]
    else:
        strain = args.input.split('.')[0]

    default_bin_size = 100

    check_call("""./cnvnator -root {}.root -tree {input}""".format(strain, input=args.input), shell=True)
    check_call("""./cnvnator -root {}.root -his {bin}""".format(strain, bin=default_bin_size), shell=True)
    check_call("""./cnvnator -root {}.root -stat {bin}""".format(strain, bin=default_bin_size), shell=True)
    check_call("""./cnvnator -root {}.root -eval {bin}""".format(strain, bin=default_bin_size), shell=True)


def load_modules():
    check_call("module load perl root/6.16.00 blast+", shell=True)

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--out_folder", required=True,
                        help="Folder to save outputs")
    parser.add_argument("-i", "--input", required=True,
                        help="Input .bam file")
    parser.add_argument("-cnv", "--cnv_path", required=True, default='/srv/scratch/lanlab/CNVnator/CNVnator_v0.3.3/src/cnvnator',
                        help="Path to cnvnator program.")


    args = parser.parse_args()

    return args

def main():
    args = parseargs()

    #load modules
    load_modules(args)

    #calculate the bin size
    cnv_bin_size = calc_cnv_bin_size(args)

    #run cnvnator

    #blast outputs - false positive checker

if __name__ == '__main__':
    main()
