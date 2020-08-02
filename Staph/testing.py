import argparse
from time import sleep as sl
import multiprocessing as mp
import pandas as pd


def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args

def process(line):

    pass

def main():
    args = parseargs()

    # # init objects
    # pool = mp.Pool(mp.cpu_count())
    # jobs = []
    #
    # # create jobs
    # with open("/Users/liamcheneyy/Desktop/out.txt") as f:
    #     for line in f:
    #         jobs.append(pool.apply_async(process(line), (line)))
    #
    # # wait for all jobs to finish
    # for job in jobs:
    #     job.get()
    #
    # # clean up
    # pool.close()



# import Bio.SeqIO
# save_list =[]
# for record in Bio.SeqIO.parse('/Users/liamcheneyy/Desktop/2020-06-18-dNdS/aligned_SACOL0001.fasta','fasta'):
# 	aa_seq = record.translate(gap='-')
# 	save_list.append(aa_seq)
# 	aa_seq.id = record.id
# 	aa_seq.description = ''
#
# Bio.SeqIO.write(save_list, '/Users/liamcheneyy/Desktop/test.fasta','fasta')
# print(save_list)
#
# from Bio.Alphabet import generic_dna, generic_protein
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from Bio.Align import MultipleSeqAlignment
# from Bio.codonalign import build
# import Bio.SeqIO
# import glob
# seqs_list = []
# for record in Bio.SeqIO.parse('/Users/liamcheney/Desktop/2020-06-18-dNdS/*.fasta','fasta'):
#     # seqs_list.append(record)


# pro_list = []
# for el in seqs_list:
#     pro = el.translate()
#     pro.id = el.id
#     pro_list.append(pro)
#     print(pro.id, len(pro.seq))
#     sl(1)
#
# pro_align = MultipleSeqAlignment(pro_list)
#
# seq1 = SeqRecord(Seq('ATGTCTCGT', alphabet=generic_dna), id='pro1')
# seq2 = SeqRecord(Seq('ATGCGT', alphabet=generic_dna), id='pro2')
# seq3 = SeqRecord(Seq('ATGTCTCGT', alphabet=generic_dna), id='pro3')
# s = [seq1,seq2,seq3]
#
# pro1 = SeqRecord(Seq('MSR', alphabet=generic_protein), id='pro1')
# pro2 = SeqRecord(Seq('M-R', alphabet=generic_protein), id='pro2')
# pro3 = SeqRecord(Seq('MSR', alphabet=generic_protein), id='pro3')
# p = [pro1,pro2,pro3]
#
# aln = MultipleSeqAlignment(p)
# codon_aln = build(aln, s)
# print(codon_aln)
# print(pro2)
# x = Bio.codonalign.codonseq.cal_dn_ds(codon_aln)
# print(x)

    # import glob
    # import Bio.SeqIO
    # save_dict = {}
    # forward_start = ['ATG','GTG','TTG']
    # reverse_sarts = ['CAT','CAC','AAT']
    #
    # tcount = 0
    # fcount = 0
    # s_list =[]
    # for filename in glob.iglob('/Users/liamcheney/Desktop/Staphylococcus/SACOL0007.fasta'):
    #
    #     for records in Bio.SeqIO.parse(filename, 'fasta'):
    #         print(records.id)
    #         print(records.seq.translate())
    #         print(str(records.seq.translate()).count('*'))
    #
    #         print(records.seq.reverse_complement().translate())
    #         print(str(records.seq.reverse_complement().translate()).count('*'))

# import glob
# from Bio import SeqIO
# infile = open('/Users/liamcheneyy/Desktop/Untitled.txt').read().splitlines()
# save_list = []
# save = []
# for filename in glob.iglob('/Users/liamcheneyy/Desktop/old_refs.fasta'):
#     for records in SeqIO.parse(filename,'fasta'):
#         locus=records.id.split(':')[0]
#         if locus in infile:
#             save.append(records)
#         if locus in infile and locus not in save_list:
#             save_list.append(locus)
#
# SeqIO.write(save, '/Users/liamcheneyy/Desktop/refs.fasta','fasta')
# print(len(save_list))
#


    # pinfile = '/Users/liamcheneyy/Downloads/assembly_result.xml'
    # infile = open(pinfile).read()
    #
    # seps = infile.split('\n\n')
    # for ind in seps[:-1]:
    #     bio = ''
    #     GCA = ''
    #     lines = ind.splitlines()
    #     for line in lines:
    #         if '<AssemblyAccession>' in line:
    #             GCA = line.split('>')[1].split('<')[0]
    #         if '<BioSampleAccn>' in line:
    #             bio = line.split('>')[1].split('<')[0]
    #
    #     print(GCA, bio)

    # import glob
    # save = []
    # for filename in glob.iglob('/Users/liamcheneyy/Downloads/genome_assemblies_wgs_gbff/ncbi-genomes-2020-07-20/*gbff'):
    #     acc= filename.split('/')[-1].split('.')[0]
    #     infile=open(filename).read().splitlines()
    #     country = ''
    #     year = ''
    #     bio = ''
    #     for line in infile:
    #         if 'country' in line:
    #             country = line.split('=')[1].strip('"').upper()
    #             if ':' in country:
    #                 country = country.split(':')[0].strip('"').upper()
    #
    #         if 'BioSample:' in line:
    #             bio = line.split(':')[1].strip().upper()
    #             print(bio)
    #
    #         if 'collection_date' in line:
    #             year = line.split('=')[1].strip('"')
    #             if '-' in year:
    #                 inds  = year.split('-')
    #                 for i in inds:
    #                     if len(i) == 4:
    #                         year = i
    #     save.append([bio, acc,year, country])
    #
    # with open('/Users/liamcheneyy/Desktop/out.txt','w') as out:
    #     for i in save:
    #         for q in i:
    #             out.write(q + '\t')
    #         out.write('\n')

#     save_list = []
    #     locus = filename.split('/')[-1].strip('.fasta')
    #     count = 0
    #     forward = False
    #     for records in Bio.SeqIO.parse(filename, 'fasta'):
    #         if '-' not in records.id:
    #             if records.seq[0:3] not in forward_start and records.seq[-3:] not in reverse_sarts:
    #                 s_list.append(records.id)
    #
    # s_list = list(set(s_list))
    # for n in s_list:
    #     print(n)

    df = pd.read_csv('/Users/liamcheneyy/Desktop/MGT_isolate_data.txt', sep='\t')

    for col in df:
        if 'MGT' in col:
            sub = df[col]
            uniques = sub.unique()
            print(uniques)




if __name__ == '__main__':
    main()
