import argparse
from time import sleep as sl
import glob
import pandas as pd

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args

def extracting_alleles(input_folder):

    save_dict = {}
    for filename in glob.iglob(input_folder + "/*alleles.fasta"):
        infile = open(filename,'r').read().splitlines()
        accession = filename.split('/')[-1].strip('_alleles.fasta').strip('.fn')
        allele_list = []
        for line in infile:
            if '>' in line and '7_gene_ST:' not in line:
                allele = line.strip('\n').strip('>')
                allele_list.append(allele)
        save_dict[accession] = allele_list

    loci_list = []

    first_strain = list(save_dict.keys())[0]
    for i in save_dict[first_strain]:
        locus = i.split(':')[0]
        loci_list.append(locus)

    return save_dict, loci_list

def calculate_missing_info(genome_alleles_dict):
    save_dict = {}

    for key, value in genome_alleles_dict.items():
        new_allele_count = 0
        worked_allele_count = 0
        exact_hit_count = 0
        missing_count = 0
        failed_filter_count = 0
        noblast_count = 0
        duplication_count = 0
        too_much_missing_count = 0
        too_long_count = 0
        for i in value:

            #handling missing alleles
            if ':0' in i:
                missing_count = missing_count + 1

                if 'failed_filter' in i:
                    failed_filter_count = failed_filter_count + 1

                elif 'possible_duplication' in i:
                    duplication_count = duplication_count + 1

                elif 'no_blast_hits' in i:
                    noblast_count = noblast_count + 1

                elif 'unscorable_too_much_missing' in i:
                    too_much_missing_count = too_much_missing_count + 1

                elif 'unscorable_too_long' in i:
                    too_long_count = too_long_count + 1

            #handling worked alleles
            if ':0' not in i:
                worked_allele_count = worked_allele_count + 1

                if 'new' in i:
                    new_allele_count = new_allele_count + 1

                elif ':1' in i:
                    exact_hit_count = exact_hit_count + 1


        save_dict[key] ={'worked':worked_allele_count, 'new':new_allele_count, 'exact':exact_hit_count, 'unscorable_too_much_missing':too_much_missing_count, 'unscorable_too_long':too_long_count, 'no_blast_hits':noblast_count, 'possible_duplication':duplication_count, 'failed_filter':failed_filter_count}

    return save_dict

def calculate_missing_alleles(loci_list, genome_alleles_dict):


    #create a dict with all loci and count numbers of zero
    locus_missing_info_dict = {}
    for locus in loci_list:
        new_allele_count = 0
        worked_allele_count = 0
        exact_hit_count = 0
        missing_count = 0
        failed_filter_count = 0
        noblast_count = 0
        duplication_count = 0
        too_much_missing_count = 0
        too_long_count = 0

        locus_missing_info_dict[locus] ={'worked':worked_allele_count, 'new':new_allele_count, 'exact':exact_hit_count, 'missing_count':missing_count, 'unscorable_too_much_missing':too_much_missing_count, 'unscorable_too_long':too_long_count, 'no_blast_hits':noblast_count, 'possible_duplication':duplication_count, 'failed_filter':failed_filter_count}

    #go over each genome and append the information to existing loci dict
    for key, value in genome_alleles_dict.items():
        for i in value:
            locus = i.split(':')[0]

            if 'failed_filter' in i:
                locus_missing_info_dict[locus]['failed_filter'] += 1

            if ':0' in i:
                locus_missing_info_dict[locus]['missing_count'] += 1

            if 'possible_duplication' in i:
                locus_missing_info_dict[locus]['possible_duplication'] += 1

            if 'no_blast_hits' in i:
                locus_missing_info_dict[locus]['no_blast_hits'] += 1

            if 'unscorable_too_much_missing' in i:
                locus_missing_info_dict[locus]['unscorable_too_much_missing'] += 1

            if 'unscorable_too_long' in i:
                locus_missing_info_dict[locus]['unscorable_too_long'] += 1

            if ':1' in i:
                locus_missing_info_dict[locus]['exact'] += 1

            if ':0' not in i:
                locus_missing_info_dict[locus]['worked'] += 1

            if 'new' in i:
                locus_missing_info_dict[locus]['new'] += 1


    return locus_missing_info_dict

def pandas_processing(locus_missing_info_dict, genome_alleles_dict, percen_out):
    df = pd.DataFrame(locus_missing_info_dict)
    df = df.T
    df = df.sort_values(by=['worked'], ascending=False)
    df = df[['worked', 'new', 'exact', 'missing_count', 'failed_filter', 'no_blast_hits', 'unscorable_too_much_missing',
             'unscorable_too_long', 'possible_duplication']]

    if percen_out:
        # number of strains
        strains_number = len(genome_alleles_dict.keys())

        df = df.apply(lambda x: round((x / strains_number * 100), 2))

    return df

def main():
    args = parseargs()

#write out as percentage or ints
    percen_out = True

    #reading in paths
    input_folder = "/Users/liamcheneyy/Desktop/90_alleles/"

    #get alleles per genome
    genome_alleles_dict,loci_list = extracting_alleles(input_folder)

    #calc number of missing alleles per genome
    missing_info_dict = calculate_missing_info(genome_alleles_dict)

    #calc the how many time an allele is missing across dataset
    locus_missing_info_dict = calculate_missing_alleles(loci_list, genome_alleles_dict)

    #extra pandas processing
    out_csv = pandas_processing(locus_missing_info_dict, genome_alleles_dict, percen_out)

    #writing out
    out_path = '/Users/liamcheneyy/Desktop/90_results.csv'
    out_csv.to_csv(out_path,sep=',')

if __name__ == '__main__':
    main()
