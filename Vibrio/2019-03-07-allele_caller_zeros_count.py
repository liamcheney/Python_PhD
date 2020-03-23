from time import sleep as sl
import glob
from Bio import SeqIO
import progressbar
import os, os.path

allele_call_input_path = "/srv/scratch/lanlab/liam/2019-04-01-reads_to_alleles/all_allele_profiles/seventh_only_alleles/*_alleles.fasta"
outfile_path = "/srv/scratch/lanlab/liam/2019-04-01-reads_to_alleles/alleles_counts/"

##count the number of zero alleles for each genome
def genome_zero_count(input_path):
    print("Caculating number of zero loci per genome")
    sl(1)
    bar = progressbar.ProgressBar(maxval=1400,
                                  widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()
    genome_count = 0
    genome_zero_count_dict = {}
    for filename in glob.iglob(input_path):
        genome_count = genome_count + 1
        bar.update(genome_count)
        genome_zero_count = 0
        accession = filename.split('/')[-1].strip('_allele.fasta')
        with open(filename) as file:
            for line in file:
                if ":0" in line:
                    genome_zero_count = genome_zero_count + 1
        genome_zero_count_dict[accession] = genome_zero_count
    return genome_zero_count_dict

##count the number of zero calls for each locus across the dataset
def locus_tag_zero_count(allele_call_input_path, number_of_strains):
    print('\n')
    print("Caculating number of zero calls per loci.")
    sl(1)
    bar = progressbar.ProgressBar(maxval=number_of_strains + 1,
                                  widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()
    genome_count = 0
    zero_per_locus_dict = {}
    for filename in glob.iglob(allele_call_input_path):
        genome_count = genome_count + 1
        bar.update(genome_count)
        accession = filename.split('/')[-1].strip('_allele.fasta')
        for record in SeqIO.parse(filename,"fasta"):
            locus_tag = record.id.split(':')[0]


            if locus_tag not in zero_per_locus_dict.keys():
                zero_per_locus_dict[locus_tag] = 0
            if ":0" in str(record.id):
                zero_per_locus_dict[locus_tag] = zero_per_locus_dict[locus_tag] + 1

    return zero_per_locus_dict

##count number of negs per locus
def locus_tag_negative(allele_call_input_path, number_of_strains):
    print('\n')
    print("Caculating number of negative calls per loci.")
    sl(1)
    bar = progressbar.ProgressBar(maxval=number_of_strains + 1,
                                  widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()
    genome_count = 0
    negative_dict = {}
    for filename in glob.iglob(allele_call_input_path):
        genome_count = genome_count + 1
        bar.update(genome_count)
        accession = filename.split('/')[-1].strip('_allele.fasta')
        for record in SeqIO.parse(filename,"fasta"):
            locus_tag = record.id.split(':')[0]
            if locus_tag not in negative_dict.keys():
                negative_dict[locus_tag] = 0
            if "N" in str(record.seq):
                negative_dict[locus_tag] = negative_dict[locus_tag] + 1
    return negative_dict

## writing out to files
def print_out_main(allele_call_input_path, outfile_path):

    # genome_zero_count_dict = genome_zero_count(allele_call_input_path)
    # with open(outfile_path + 'genomes_zeros_count.csv', 'w') as outfile:
    #     for key, value in genome_zero_count_dict.items():
    #         outfile.write(key + ',' + str(value) + '\n')

    fasta_path = allele_call_input_path.split('*')[0]
    number_of_strains = len(os.listdir(fasta_path))

    zero_per_locus_dict = locus_tag_zero_count(allele_call_input_path, number_of_strains)
    with open(outfile_path + 'locus_zero_counts.csv', 'w') as outfile:
        for key, value in zero_per_locus_dict.items():
            outfile.write(key + ',' + str(value) + '\n')
    # #
    negative_dict = locus_tag_negative(allele_call_input_path, number_of_strains)
    with open(outfile_path + 'locus_negative_counts.csv', 'w') as outfile:
        for key, value in negative_dict.items():
            outfile.write(key + ',' + str(value) + '\n')

print_out_main(allele_call_input_path,outfile_path)
