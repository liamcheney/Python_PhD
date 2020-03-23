import argparse
from time import sleep as sl


def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args

def cal_stats(infile):

    strain_info = infile.split('Parsing inputs')

    for element in strain_info:
        new = element.split('Writing outputs')[-1].split('\n')
        strain = ''
        called = ''
        intact = ''
        partial = ''
        unscorable_too_much_missing = ''
        no_blast_hits = ''
        missing = ''
        duplicate = 999

        for line in new:
                if 'strain:' in line:
                    strain = line.split(':')[-1].strip(' ')

                if 'called:' in line:
                    called = int(line.split(':')[-1])

                if 'intact:' in line:
                    intact = int(line.split(':')[-1])

                if 'partial:' in line:
                    partial = int(line.split(':')[-1])

                if 'unscorable_too_much_missing:' in line:
                    unscorable_too_much_missing = int(line.split(':')[-1])

                if 'no_blast_hits:' in line:
                    no_blast_hits = int(line.split(':')[-1])

                if 'duplicate' in line:
                    no_blast_hits = int(line.split(':')[-1])

        missing = unscorable_too_much_missing + no_blast_hits

        print(strain, called, intact, partial, unscorable_too_much_missing, no_blast_hits, missing, duplicate)

def main():
    args = parseargs()

    path = '/Users/liamcheneyy/Desktop/poor_assembles_test/hsp70_snp16_blast70/output.txt'
    infile = open(path,'r').read()

    print(path.split('/')[-2])
    cal_stats(infile)

if __name__ == '__main__':
    main()
