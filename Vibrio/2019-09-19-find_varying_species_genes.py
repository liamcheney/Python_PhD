import argparse
from time import sleep as sl

def calculate_alleles(single_genes, all_species_alleles):

    #get list of all alleles
    alleles_list = []
    for line in all_species_alleles:
        if '>' in line:
            locus = line.strip('>').split(':')[0]
            alleles_list.append(locus)


    #just count occurence
    save_dict = {}
    for line in single_genes:
        count = 0
        line = line.strip().strip('\n')
        for element in alleles_list:
            if element == line:
                count = count + 1
        save_dict[line] = count

    return save_dict
def out_write(count_dict, outfile_path):

    with open(outfile_path,'w') as out:
        out.write('Locus Tag' + '\t' + 'Species Allele Count' + '\n')
        for key, value in count_dict.items():
            out.write(key + '\t' + str(value) + '\n')


def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args


def main():
    args = parseargs()

    single_genes_path = "/Users/liamcheneyy/Desktop/single_allele_loci.txt"
    single_genes = open(single_genes_path,'r').read().splitlines()

    all_species_alleles_path = "/Users/liamcheneyy/Desktop/all_ref_alleles.fasta"
    all_species_alleles = open(all_species_alleles_path,'r').read().splitlines()

    outfile_path = '/Users/liamcheneyy/Desktop/species_alleles_counts.txt'

    count_dict = calculate_alleles(single_genes, all_species_alleles)

    out_write(count_dict, outfile_path)
if __name__ == '__main__':
    main()
