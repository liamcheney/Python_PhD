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

    # trf_in = open('/Users/liamcheneyy/Desktop/Book5.txt').read().splitlines()
    # trf_list = []
    # for line in trf_in:
    #     col = line.split('--')
    #     s = int(col[0])
    #     e = int(col[1])
    #     ls = [s,e]
    #     trf_list.append(ls)
    #
    # meta = open('/Users/liamcheneyy/Desktop/gene_metda.txt').read().splitlines()
    # for line in meta[1:]:
    #     TR = False
    #     col = line.split('\t')
    #     locus = col[0]
    #     s = int(col[1])
    #     e = int(col[2])
    #
    #     for el in trf_list:
    #         trf_s = el[0]
    #         trf_e = el[1]
    #
    #         if (trf_s >= s) and (trf_e <= e):
    #             # print(locus, trf_s,trf_e, s, e)
    #             TR = True
    #
    #
    #     print(locus, TR, sep='\t')


if __name__ == '__main__':
    main()
