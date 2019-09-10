import argparse
import pandas as pd
import glob
from time import sleep as sl

def create_alleles_dict(mgt9_alleles):
    # print('Reading in alleles for ' + str(len(mgt9_alleles)-1) + ' strains.')

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
            if '-' not in cell:
                cell_locus = locus + '_' + str(cell)
                alleles_list.append(cell_locus)

            # elif '-' in cell:
            #     cell_locus = locus + '_' + str(cell)
            #     alleles_list.append(cell_locus)

            alleles_dict[strain] = alleles_list

    return alleles_dict

def count_alleles_per_loci(mgt9_alleles_path):

    df = pd.read_csv(mgt9_alleles_path, sep='\t', index_col=False)

    allele_count = {}
    for column in df:
        if '#' not in column:
            #not using
            # unique_values = list(df[column].unique())
            # alleles_num = len(unique_values)
            col_freqs = df[column].value_counts().to_dict()
            allele_count[column] = col_freqs
            # print(column, allele_count[column])
            # sl(1)

    return allele_count

def read_in_alleles_per_loci(alleles_in):

    alleles_dict = {}
    for line in alleles_in:
        col = line.split('\t')
        alleles_dict[col[0]] = col[1]

    return alleles_dict

def temp_write_out(alleles_dict):
    with open('/Users/liamcheneyy/Desktop/temp_file.csv','w') as out:
        for key, value in alleles_dict.items():
            out.write(key + ',')

            for i in alleles_dict[key]:
                out.write(i + '_' + alleles_dict[key][i] + ',')

            out.write('\n')

def terminating_genes_per_clade(input_genomes, strains_close_to_clade, alleles_dict, calc_alleles_strains_per_loc):

    #MAIN: calculate the alleles specific to the input clade

    #get all alleles for a clade
    shared_alleles_from_input_clade = alleles_per_clade(alleles_dict, input_genomes)

    #find the alleles in non input strains
    allele_in_all_list = alleles_from_all_strains(alleles_dict, strains_close_to_clade, input_genomes)

    #find alleles only changing in certain clade
    clade_specific_genes = alleles_specific_to_clade(allele_in_all_list, shared_alleles_from_input_clade, calc_alleles_strains_per_loc)

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
    # print("Finding Shared Alleles.")

    # print(len(alleles_dict.keys()))

    #remove input genomes and others in clade from alleles_dict
    all_remove_list = list(set(input_genomes + strains_close_to_clade))
    for element in all_remove_list:
        del alleles_dict[element]

    # print(len(alleles_dict.keys()))

    #list of all alleles for a clade
    all_alleles = []

    for key in alleles_dict.keys():
        for cell in alleles_dict[key]:
            all_alleles.append(cell)

    all_alleles = list(set(all_alleles))

    return all_alleles

def alleles_specific_to_clade(allele_in_all_list, shared_alleles_from_input_clade, calc_alleles_strains_per_loc):
    # print('Finding Clade Specific Genes.')

    #find specific genes
    specific_genes = []
    for allele in shared_alleles_from_input_clade:
        if allele not in allele_in_all_list:
            specific_genes.append(allele)
            # print(allele)

    # print(str(len(specific_genes)) + " genes Specific to input clade.")
    for i in specific_genes:
        locus = i.split('_')[0]
        alleles_num = len(calc_alleles_strains_per_loc[locus].keys())
        print(str(i) + '\t' + str(alleles_num) + '\t' + str(calc_alleles_strains_per_loc[locus]))

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args

def main():
    args = parseargs()
    input_genomes_path = '/Users/liamcheneyy/Desktop/input_clade.txt'
    input_genomes = open(input_genomes_path, 'r').read().splitlines()

    mgt9_alleles_path = '/Users/liamcheneyy/Desktop/MGT9_allele_profiles.tsv'
    mgt9_alleles = open(mgt9_alleles_path,'r').read().splitlines()

    all_mgt9_alleles_path = mgt9_alleles_path

    all_of_interest_path = '/Users/liamcheneyy/Desktop/remove_clade.txt'
    strains_close_to_clade = open(all_of_interest_path, 'r').read().splitlines()

    alleles_nums_in = '/Users/liamcheneyy/Desktop/alleles_per_loci.txt'
    alleles_in = open(alleles_nums_in, 'r').read().splitlines()

    #turn input into alleles dictionary
    #dict key:strain, value:alleles
    alleles_dict = create_alleles_dict(mgt9_alleles)

    #calculate alleles per loci, and percentage of strains each allele holds
    calc_alleles_strains_per_loc = count_alleles_per_loci(all_mgt9_alleles_path)

    #find genes for a certain TERMINATING clade
    input_clade_specific_genes = terminating_genes_per_clade(input_genomes, strains_close_to_clade, alleles_dict, calc_alleles_strains_per_loc)


if __name__ == '__main__':
    main()
