from time import sleep as sl
import glob
from Bio import SeqIO

reference_cds_in = "/Users/liamcheneyy/Desktop/ref/GCA_000006745.1_ASM674v1_cds_from_genomic.fna"
mgt_schemes_path = "/Users/liamcheneyy/Desktop/Schemes/"
out_locus_gene_sizes_path = '/Users/liamcheneyy/Desktop/mgt_schemes_genes.csv'

in_locus_ds_path = '/Users/liamcheneyy/Desktop/ds_all_loci.txt'
out_locus_ds_path = '/Users/liamcheneyy/Desktop/mgt_schemes_ds.csv'

in_locus_dnds_path = "/Users/liamcheneyy/Desktop/dnds_all_loci.txt"
out_locus_dnds_path = '/Users/liamcheneyy/Desktop/mgt_schemes_dnds.csv'


def create_cds_dict(reference_cds_in):
    cds_dict = {}
    with open(reference_cds_in) as file:
        for record in SeqIO.parse(file,'fasta'):
            locus_tag = record.description.split()[1].split('=')[-1].strip(']').split('_')[0] + record.description.split()[1].split('=')[-1].strip(']').split('_')[1]
            cds_dict[locus_tag] = str(record.seq)
    return cds_dict

def check_gene_size(reference_cds_in, mgt_schemes_path, outfile):
    cds_dict = create_cds_dict(reference_cds_in)
    with open(outfile, 'w') as out:
        for scheme in glob.iglob(mgt_schemes_path + '*'):
            scheme_name = scheme.split('/')[-1].split('_')[0]
            out.write(scheme_name + '\n')
            infile = open(scheme, 'r').read().splitlines()
            for line in infile:
                if line in cds_dict.keys():
                    gene_length = len(cds_dict[line])
                    out.write(line + ',' + str(gene_length) + '\n')
                else:
                    out.write(line + ',' + 'No Length' + '\n')
            out.write('\n')
    out.close()

def check_ds_per_scheme(mgt_schemes_path, locus_ds_path, outfile):

    locus_ds_dict = {}
    locus_ds_list = list([x for x in open(locus_ds_path).read().splitlines()])
    for i in locus_ds_list:
        col = i.split()
        locus_ds_dict[col[0]] = col[1]

    with open(outfile, 'w') as out:
        for scheme in glob.iglob(mgt_schemes_path + '*'):
            scheme_name = scheme.split('/')[-1].split('_')[0]
            out.write(scheme_name + '\n')

            infile = open(scheme, 'r').read().splitlines()
            for line in infile:
                if line in locus_ds_dict.keys():
                    out.write(line +','+ str(locus_ds_dict[line]) + '\n')
            out.write('\n')

def check_dnds_per_scheme(mgt_schemes_path, locus_ds_path, outfile):

    locus_ds_dict = {}
    locus_ds_list = list([x for x in open(locus_ds_path).read().splitlines()])
    for i in locus_ds_list:
        col = i.split()
        locus_ds_dict[col[0]] = col[1]

    with open(outfile, 'w') as out:
        for scheme in glob.iglob(mgt_schemes_path + '*'):
            scheme_name = scheme.split('/')[-1].split('_')[0]
            out.write(scheme_name + '\n')

            infile = open(scheme, 'r').read().splitlines()
            for line in infile:
                if line in locus_ds_dict.keys():
                    out.write(line +','+ str(locus_ds_dict[line]) + '\n')
            out.write('\n')


check_ds_per_scheme(mgt_schemes_path, in_locus_ds_path, out_locus_ds_path)

check_dnds_per_scheme(mgt_schemes_path, in_locus_dnds_path, out_locus_dnds_path)

check_gene_size(reference_cds_in, mgt_schemes_path, out_locus_gene_sizes_path)

