import argparse
import pandas as pd

from time import sleep as sl
import glob

def create_alleles_dict(mgt9_alleles, input_genomes):
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

    # #check if any input strains are not in represntative MGT9 alleles
    # input_all_mgt9_path = '/Users/liamcheneyy/Desktop/vcseventh_15/grapetree/all_MGT9_allele_profiles.tsv'
    # all_mgt9_alleles = open(input_all_mgt9_path,'r').read().splitlines()
    # for el in input_genomes:
    #     if el not in alleles_dict.keys():
    #         print(el)
    #         locus_line = mgt9_alleles[0].split('\t')
    #         for line in mgt9_alleles[1:]:
    #             col = line.split('\t')
    #             print(col)
    #             sl(1)
    #             if col[0] == el:
    #                 print(col)
    #                 col_count = 3
    #                 beggining_col = 3
    #                 strain = col[0]
    #                 alleles_list = []
    #
    #                 for cell in col[beggining_col:]:
    #                     locus = locus_line[col_count]
    #                     col_count = col_count + 1
    #                     if '-' not in cell:
    #                         cell_locus = locus + '_' + str(cell)
    #                         alleles_list.append(cell_locus)
    #
    #                     alleles_dict[strain] = alleles_list
    #                     print(strain, alleles_list)

    return alleles_dict

def add_extra_alleles_profiles(mgt9_alleles, input_genomes, all_mgt9_alleles):

    strains_from_allels_in = []
    for line in mgt9_alleles[1:]:
        col = line.split('\t')
        strains_from_allels_in.append(col[0])

    #need to get list
    get_list = []
    for line in input_genomes:
        if line not in strains_from_allels_in:
            get_list.append(line)

    for line in all_mgt9_alleles[1:]:
        col = line.split('\t')
        if col[0] in get_list:
            mgt9_alleles.append(line)

    return mgt9_alleles

def count_alleles_per_loci(mgt9_alleles_path, input_clade_size):

    df = pd.read_csv(mgt9_alleles_path, sep='\t', index_col=False)

    allele_count = {}
    strains_num = df.shape[0]
    for column in df:
        if '#' not in column:
            #not using
            # unique_values = list(df[column].unique())
            # alleles_num = len(unique_values)
            # print(column, unique_values, alleles_num)

            col_freqs = df[column].value_counts().to_dict()
            allele_count[column] = col_freqs

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

def terminating_genes_per_clade(input_genomes, alleles_dict, calc_alleles_strains_per_loc, vibrio_chol_genes):

    #MAIN: calculate the alleles specific to the input clade

    #get all alleles for a clade
    shared_alleles_from_input_clade = alleles_per_clade(alleles_dict, input_genomes)

    #find the alleles in non input strains
    allele_in_all_list = alleles_from_all_strains(alleles_dict, input_genomes)

    #find alleles only changing in certain clade
    clade_specific_genes = alleles_specific_to_clade(allele_in_all_list, shared_alleles_from_input_clade, calc_alleles_strains_per_loc, vibrio_chol_genes)

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

def alleles_from_all_strains(alleles_dict, input_genomes):
    # print("Finding Shared Alleles.")

    # print(len(alleles_dict.keys()))

    #remove input genomes and others in clade from alleles_dict
    all_remove_list = list(set(input_genomes))
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

def alleles_specific_to_clade(allele_in_all_list, shared_alleles_from_input_clade, calc_alleles_strains_per_loc, vibrio_chol_genes):
    # print('Finding Clade Specific Genes.')

    #find specific genes
    specific_genes = []
    for allele in shared_alleles_from_input_clade:
        if allele not in allele_in_all_list:
            specific_genes.append(allele)
            # print(allele)

    print(str(len(specific_genes)) + " genes Specific to input clade.")
    out_list = []
    for i in specific_genes:
        locus = i.split('_')[0]
        alleles_num = len(calc_alleles_strains_per_loc[locus].keys())

        check_list = False
        if locus in vibrio_chol_genes:
            check_list = True

        out_list.append([str(i), str(alleles_num), str(calc_alleles_strains_per_loc[locus]), str(check_list)])
        # print(str(i) + '\t' + str(alleles_num) + '\t' + str(calc_alleles_strains_per_loc[locus]))

    #sort list by number of alleles
    out_list.sort(key=lambda x:(int(x[1]), x[-1]), reverse=False)
    for el in out_list:
        out_str = '\t'.join(el)
        print(out_str)

    print()

def work_flow(input_genomes_path, mgt9_alleles_path, all_mgt9_alleles_path, vibrio_core_list_path):

    input_genomes = open(input_genomes_path, 'r').read().splitlines()

    mgt9_alleles = open(mgt9_alleles_path, 'r').read().splitlines()

    all_mgt9_alleles = open(all_mgt9_alleles_path, 'r').read().splitlines()

    vibrio_chol_genes = open(vibrio_core_list_path, 'r').read().splitlines()

    # if comparing strains not in original allele alignment
    mgt9_alleles = add_extra_alleles_profiles(mgt9_alleles, input_genomes, all_mgt9_alleles)

    # turn input into alleles dictionary
    # dict key:strain, value:alleles
    alleles_dict = create_alleles_dict(mgt9_alleles, input_genomes)

    # calculate alleles per loci, and percentage of strains each allele holds
    input_clade_size = len(input_genomes)
    calc_alleles_strains_per_loc = count_alleles_per_loci(all_mgt9_alleles_path, input_clade_size)

    # find genes for a certain TERMINATING clade
    input_clade_specific_genes = terminating_genes_per_clade(input_genomes, alleles_dict, calc_alleles_strains_per_loc,
                                                             vibrio_chol_genes)


def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args

def main():
    args = parseargs()

    input_genomes_path = '/Users/liamcheneyy/Desktop/input_clades/'
    mgt9_alleles_path = '/Users/liamcheneyy/Desktop/MGT9_allele_profiles.tsv'
    all_mgt9_alleles_path = '/Users/liamcheneyy/Desktop/vcseventh_21/grapetree/seventh/all_MGT8_allele_profiles.tsv'
    vibrio_core_list_path = '/Users/liamcheneyy/Desktop/MGT8_gene_accessions.txt'

    loop = True
    if loop:
        for filename in glob.iglob(input_genomes_path + '/*'):
            file_name = filename.split('/')[-1].split('.')[0]
            print(file_name)
            work_flow(filename, mgt9_alleles_path, all_mgt9_alleles_path, vibrio_core_list_path)
    else:
        work_flow(input_genomes_path, mgt9_alleles_path, all_mgt9_alleles_path, vibrio_core_list_path)


if __name__ == '__main__':
    main()
