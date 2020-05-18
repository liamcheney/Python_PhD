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


    # ###go over allele profiles and get summary stats
    # genome_zero_neg_dict = {}
    # gene_zero_dict = {}
    # for filename in glob.iglob('/Users/liamcheneyy/Desktop/alleles/*fasta'):
    #     bios = SeqIO.parse(filename,'fasta')
    #     accession = filename.split('/')[-1].split('_')[0]
    #     print(accession)
    #
    #     genome_zero_neg_dict[accession] = {'zero':0,'neg':0}
    #
    #     ##count zeros per genome
    #     for record in bios:
    #         gene = record.id.split(':')[0]
    #
    #         if gene not in gene_zero_dict.keys():
    #             gene_zero_dict[gene] = {'zero':0, 'intact':0, 'new':0, 'neg':0}
    #
    #         if ":0" in record.id:
    #             gene_zero_dict[gene]['zero'] = gene_zero_dict[gene]['zero'] + 1
    #             genome_zero_neg_dict[accession]['zero'] = genome_zero_neg_dict[accession]['zero'] + 1
    #
    #         if ":1" in record.id or ":2" in record.id or ":3" in record.id or ":4" in record.id or ":5" in record.id or ":6" in record.id or ":7" in record.id:
    #             gene_zero_dict[gene]['intact'] = gene_zero_dict[gene]['intact'] + 1
    #
    #         if ":new" in record.id:
    #             gene_zero_dict[gene]['new'] = gene_zero_dict[gene]['new'] + 1
    #
    #         if 'N' in record.seq:
    #             gene_zero_dict[gene]['neg'] = gene_zero_dict[gene]['neg'] + 1
    #             genome_zero_neg_dict[accession]['neg'] = genome_zero_neg_dict[accession]['neg'] + 1
    #
    #
    # with open('/Users/liamcheneyy/Desktop/alleles/genome_zeros.txt', 'w') as out:
    #     out.write('Accessions' + '\t' + 'Zero' + '\t' + 'Neg' + '\n')
    #     for key, value in genome_zero_neg_dict.items():
    #         out.write(key + '\t')
    #         for i in value:
    #             out.write(str(value[i]) + '\t')
    #         out.write('\n')
    #
    # with open('/Users/liamcheneyy/Desktop/alleles/gene_zeros.txt', 'w') as out:
    #     out.write('Gene' + '\t' + 'Zero' + '\t' + 'Intact' + '\t' + 'New' + '\t' + 'Neg' + '\n')
    #     for key, value in gene_zero_dict.items():
    #         out.write(key + '\t')
    #         for i in value:
    #             out.write(str(value[i]) + '\t')
    #         out.write('\n')




        # genome_zero_dict[accession] = genome_zero_conut

if __name__ == '__main__':
    main()