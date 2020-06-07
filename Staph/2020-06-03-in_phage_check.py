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

    regions = [[354767,398136], [900796,922004]]

    infile = open('/Users/liamcheneyy/Desktop/Book4.txt','r').read().splitlines()

    for line in infile:
        phage = False
        col = line.split('\t')
        s = int(col[1])
        e = int(col[2])
        locus = col[0]

        for el in regions:
            el_s = el[0]
            el_e = el[1]

            if (s >= el_s) and (e <= el_e):
                print(locus, el, s, e)


if __name__ == '__main__':
    main()
