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

    recom_snps = open('/Users/liamcheneyy/Desktop/2020-06-18-recombination/norecom.rec.location').read().splitlines()

    infile = open('/Users/liamcheneyy/Desktop/2020-06-18-recombination/lociLocationsInRef.txt').read().splitlines()

    save_list = []
    for SNP in recom_snps:
        # print(SNP)
        iSNP = int(SNP)
        for line in infile:
            col = line.split('\t')
            l = col[0]
            s = int(col[1])
            e = int(col[2])

            if (iSNP >= s) and (iSNP <= e):
                save_list.append(l)

    save_list = list(set(save_list))

    for el in save_list:
        print(el)

if __name__ == '__main__':
    main()
