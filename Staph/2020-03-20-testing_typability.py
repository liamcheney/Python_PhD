import argparse
from time import sleep as sl
from Bio import SeqIO
import glob

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args


def main():
    args = parseargs()

    #go over allele profiles and get summary stats
    genome_zero_neg_dict = {}
    gene_zero_dict = {}
    for filename in glob.iglob('/Users/liamcheney/Desktop/alleles/*fasta'):
        bios = SeqIO.parse(filename,'fasta')
        accession = filename.split('/')[-1].split('_')[0]
        print(accession)

        genome_zero_neg_dict[accession] = {'zero':0,'neg':0}

        ##count zeros per genome
        for record in bios:
            gene = record.id.split(':')[0]

            if gene not in gene_zero_dict.keys():
                gene_zero_dict[gene] = {'zero':0, 'intact':0, 'new':0, 'neg':0}

            if ":0" in record.id:
                gene_zero_dict[gene]['zero'] = gene_zero_dict[gene]['zero'] + 1
                genome_zero_neg_dict[accession]['zero'] = genome_zero_neg_dict[accession]['zero'] + 1

            if record.id.split(':')[-1].isdigit():
                gene_zero_dict[gene]['intact'] = gene_zero_dict[gene]['intact'] + 1

            if ":new" in record.id:
                gene_zero_dict[gene]['new'] = gene_zero_dict[gene]['new'] + 1

            if 'N' in record.seq:
                gene_zero_dict[gene]['neg'] = gene_zero_dict[gene]['neg'] + 1
                genome_zero_neg_dict[accession]['neg'] = genome_zero_neg_dict[accession]['neg'] + 1


    with open('/Users/liamcheney/Desktop/alleles/genome_zeros.txt', 'w') as out:
        out.write('Accessions' + '\t' + 'Zero' + '\t' + 'Neg' + '\n')
        for key, value in genome_zero_neg_dict.items():
            out.write(key + '\t')
            for i in value:
                out.write(str(value[i]) + '\t')
            out.write('\n')

    with open('/Users/liamcheney/Desktop/alleles/gene_zeros.txt', 'w') as out:
        out.write('Gene' + '\t' + 'Zero' + '\t' + 'Intact' + '\t' + 'New' + '\t' + 'Neg' + '\n')
        for key, value in gene_zero_dict.items():
            out.write(key + '\t')
            for i in value:
                out.write(str(value[i]) + '\t')
            out.write('\n')

    # ##remove locus from alleles based on input list
    # infile = open('/Users/liamcheneyy/Desktop/xx.txt','r').read().splitlines()
    # save_dict  = {}
    # for filename in glob.iglob('/Users/liamcheneyy/Desktop/fixed_alleles/*fasta'):
    #     records = SeqIO.parse(filename, 'fasta')
    #
    #     accession = filename.split('/')[-1]
    #     save_dict[accession] = {}
    #     for record in records:
    #         if not record.id.split(':')[0].strip('>') in infile:
    #             save_dict[accession][record.id] = record.seq
    #
    # for k, v in save_dict.items():
    #     print(k)
    #     with open('/Users/liamcheneyy/Desktop/outs/' + k ,'w') as out:
    #         for id,seq in v.items():
    #             out.write(">" + id + "\n")
    #
    #             if "new" in id:
    #                 out.write(str(seq) + "\n\n")
    #             elif not "new" in id:
    #                 out.write(str(seq) + "\n")


    # # remove locus from reference file
    # keep = open('/Users/liamcheneyy/Desktop/out.txt').read().splitlines()
    #
    # records = SeqIO.parse('/Users/liamcheneyy/Desktop/refs.fasta','fasta')
    #
    # rem_list = []
    # save_list = []
    # for record in records:
    #     locus = record.id.split(':')[0]
    #
    #     if locus in keep:
    #         save_list.append(record)
    #         rem_list.append(locus)


    # SeqIO.write(save_list,'/Users/liamcheneyy/Desktop/fixrefs.fasta','fasta')
    # for i in save_list:
    #     print(i)

if __name__ == '__main__':
    main()