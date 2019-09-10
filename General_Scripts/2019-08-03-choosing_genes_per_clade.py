import argparse
from time import sleep as sl

def create_alleles_dict(mgt9_alleles):
    print('Reading in alleles for ' + str(len(mgt9_alleles)-1) + ' strains.')

    alleles_dict = {}
    locus_line = mgt9_alleles[0].split('\t')
    for line in mgt9_alleles[1:]:
        col = line.split('\t')
        col_count = 3
        beggining_col = 3
        strain = col[0]
        alleles_list = []

        for cell in col[beggining_col:]:
            locus = locus_line[col_count]
            col_count = col_count + 1
            if '-' not in cell and int(cell) > 1:
                cell_locus = locus + '_' + str(cell)
                alleles_list.append(cell_locus)

            elif '-' in cell:
                cell_locus = locus + '_' + str(cell)
                alleles_list.append(cell_locus)

            alleles_dict[strain] = alleles_list

    return alleles_dict

def temp_write_out(alleles_dict):
    with open('/Users/liamcheneyy/Desktop/temp_file.csv','w') as out:
        for key, value in alleles_dict.items():
            out.write(key + ',')

            for i in alleles_dict[key]:
                out.write(i + '_' + alleles_dict[key][i] + ',')

            out.write('\n')

def terminating_genes_per_clade(input_genomes, strains_close_to_clade, alleles_dict):

    #then which of those genes not found in another strains

    #get all alleles for a clade
    shared_alleles_from_input_clade = alleles_per_clade(alleles_dict, input_genomes)

    #find the alleles in non input strains
    allele_in_all_list = alleles_from_all_strains(alleles_dict, strains_close_to_clade, input_genomes)

    #find alleles only changing in certain clade
    clade_specific_genes = alleles_specific_to_clade(allele_in_all_list, shared_alleles_from_input_clade)

    return clade_specific_genes

def alleles_per_clade(alleles_dict, input_genomes):

    #make subset for input strains
    sub_dict = {}
    for key in alleles_dict.keys():
        if key in input_genomes:
            sub_dict[key] = alleles_dict[key]

    #find the genes which are common to all strains
    #make list of all genes
    all_genes_lists = []
    for key in sub_dict.keys():
        all_genes_lists.append(sub_dict[key])

    #make list of commonly shared genes
    commonly_shared_genes = list(set.intersection(*map(set,all_genes_lists)))

    return commonly_shared_genes

def alleles_from_all_strains(alleles_dict, input_genomes, strains_close_to_clade):
    print("Finding Shared Alleles.")

    #remove input genomes and others in clade from alleles_dict
    all_remove_list = list(set(input_genomes + strains_close_to_clade))
    for element in all_remove_list:
        print(element)
        del alleles_dict[element]


    #list of all alleles for a clade
    all_alleles = []

    for key in alleles_dict.keys():
        for cell in alleles_dict[key]:
            all_alleles.append(cell)

    all_alleles = list(set(all_alleles))

    return all_alleles

def alleles_specific_to_clade(allele_in_all_list, shared_alleles_from_input_clade):
    print('Finding Clade Specific Genes.')

    #find specific genes
    specific_genes = []
    for allele in shared_alleles_from_input_clade:
        if allele not in allele_in_all_list:
            specific_genes.append(allele)
            # print(allele)

    print(str(len(specific_genes)) + " genes Specific to input clade.")
    for i in specific_genes:
        print(i)

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args

def main():
    args = parseargs()

    # input_genomes_path = '/Users/liamcheneyy/Desktop/small_input.txt'
    # mgt9_alleles_path = '/Users/liamcheneyy/Desktop/small_alleles.tsv'

    input_genomes_path = '/Users/liamcheneyy/Desktop/small.txt'
    input_genomes = open(input_genomes_path,'r').read().splitlines()

    mgt9_alleles_path = '/Users/liamcheneyy/Desktop/MGT9_allele_profiles.tsv'
    mgt9_alleles = open(mgt9_alleles_path,'r').read().splitlines()

    all_of_interest_path = '/Users/liamcheneyy/Desktop/remove_clade.txt'
    strains_close_to_clade = open(all_of_interest_path, 'r').read().splitlines()

    #dict key:strain, value:alleles
    alleles_dict = create_alleles_dict(mgt9_alleles)

    #find genes for a certain TERMINATING clade
    input_clade_specific_genes = terminating_genes_per_clade(input_genomes, strains_close_to_clade, alleles_dict)


if __name__ == '__main__':
    main()
