import argparse
from time import sleep as sl
import glob

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args

def get_info(infile, wanted):

    results = []
    for line in infile:
        if (wanted in line) and ('serogroup_cholerae' not in line) and ('biotype_cholerae' not in line):
            blast_coverage = float(line.split('\t')[2])
            gene = line.split('\t')[1]
            results.append([gene, blast_coverage])

    results.sort(key=lambda x: x[1], reverse=True)

    if len(results) > 0:
        # top_hit = results[0][0]
        return True
    else:
        return False


def main():
    args = parseargs()
    input_path = '/Users/liamcheneyy/Desktop/cholera_finder'
    out_path = '/Users/liamcheneyy/Desktop/'

    ###creating list for all genes
    wanted_db = []
    for filename in glob.iglob(input_path + '/*/results_tab.tsv'):
        infile = open(filename,'r').read().splitlines()
        for line in infile[1:]:
            col = line.split('\t')
            gene = col[1]
            if ('serogroup_cholerae' not in line) and ('biotype_cholerae' not in line):
                if gene not in wanted_db:
                    wanted_db.append(gene)

    ###collecting metadata from all strains
    all_results = {}
    all_strains = []
    for filename in glob.iglob(input_path + '/*/results_tab.tsv'):
        infile = open(filename,'r').read().splitlines()
        strain = filename.split('/')[-2]
        all_strains.append(strain)

        all_results[strain] = {}
        for element in wanted_db:
            all_results[strain][element] = ''

        for wanted in wanted_db:
            all_results[strain][wanted] = get_info(infile, wanted)

    ###writing out metadata
    with open(out_path + '/choleraefinder_results.csv', 'w') as out:
        out.write('Strain' + ',')
        for element in wanted_db:
            out.write(element + ',')
        out.write('\n')

        for strain in all_strains:
            out.write(strain + ',')
            for element in wanted_db:
                out.write(str(all_results[strain][element]) + ',')
            out.write('\n')

if __name__ == '__main__':
    main()
