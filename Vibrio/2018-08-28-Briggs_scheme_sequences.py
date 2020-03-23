from time import sleep as sl
import glob
from Bio import SeqIO
import os

infile = '/Users/liam/Desktop/Vibrio_D3A_loci_accession_lists/'
ref_alleles = '/Users/liam/Desktop/Vibrio_D3A_alleles_ref.fasta'
output = '/Users/liam/Desktop/individ_alleles/'

def brigs_scheme_to_sequence(MGT_schemes_folder, ref_alleles_fasta, output):

    ref_alleles_dict = {}
    ref_alleles_in = SeqIO.parse(ref_alleles_fasta, 'fasta')
    for allele in ref_alleles_in:
        ref_alleles_dict[allele.id[:-2]] = allele.seq

    for filename in glob.iglob(MGT_schemes_folder + "*.txt"):
        scheme = filename.split('/')[-1].split('_')[0]
        os.mkdir(output + scheme)

        scheme_list = []
        with open(filename, 'r') as file:
            for line in file:
                scheme_list.append(line.rstrip('\n'))

        for element in scheme_list:
            # with open(output + scheme + '.fasta', 'a') as outfile:
            with open(output + scheme + '/' + element + '.fasta', 'a') as outfile:
                outfile.write(">" + element + '\n')
                outfile.write(str(ref_alleles_dict[element]) + '\n')

brigs_scheme_to_sequence(infile, ref_alleles, output)

