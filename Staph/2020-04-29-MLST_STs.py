import argparse
from time import sleep as sl


def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args


def main():
    args = parseargs()

    ##create a local database of novel MLST alleles and select a strain from each

    infile = open('/Users/liamcheneyy/Desktop/Book2.txt').read().splitlines()

    combs = []
    combs_dict = {}
    count = 1
    for line in infile:
        col = line.split('\t')
        alleles = col[3:]
        comb_str = '-'.join(alleles)

        if comb_str not in combs_dict.keys():
            combs_dict[comb_str] = "n" + str(count)
            print(line, combs_dict[comb_str])
            count = count + 1

        elif comb_str in combs_dict.keys():
            print(line, combs_dict[comb_str])



if __name__ == '__main__':
    main()
