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

if __name__ == '__main__':
    main()


#figure what the ctxb alleles are


listn = [1,2,3,4,5,6,7,1,2,3]

print(listn)
