import argparse
from time import sleep as sl
import pandas

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args


def main():
    args = parseargs()

    infile_path = '/Users/liamcheneyy/Desktop/RecDetect/vc_species.real.withN.tsv'
    infile = open(infile_path,'r').read().splitlines()

    strains_num = len(infile[0].split('\t')) - 1
    save_dict = {}

    for i in range(1, strains_num + 1, 1):
        strain = infile[0].split('\t')[i]
        print(strain)
        save_list = []
        for line in infile[1:]:
            col = line.split('\t')[i]
            save_list.append(col)

        complete = ''.join(save_list)
        save_dict[strain] = complete

    out_path = '/Users/liamcheneyy/Desktop/RecDetect/vc_species.real.withN.fasta'
    with open(out_path,'w') as out:
        for k, v in save_dict.items():
            out.write('>' + k + '\n')
            out.write(v + '\n')









if __name__ == '__main__':
    main()
