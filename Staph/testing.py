import argparse
from time import sleep as sl
import multiprocessing as mp


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
    # import glob
    # import numpy as np
    #
    # for filename in glob.iglob('/Users/liamcheneyy/MGT_28/Mgt/Mgt/Alleles/Staphylococcus/SACOL0001.fasta'):
    #     save_list = []
    #     diff_lengths = False
    #     for records in Bio.SeqIO.parse(filename,'fasta'):
    #         length = len(records.seq)
    #         save_list.append(length)
    #
    #         if len(records.seq) == 1361:
    #             print(records.id)
    #
    #     save_list = list(set(save_list))
    #     print(len(save_list), save_list, filename)
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

from Bio.Alphabet import generic_dna, generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.codonalign import build
import Bio.SeqIO

seqs_list = []
for record in Bio.SeqIO.parse('/Users/liamcheneyy/Desktop/2020-06-18-dNdS/SACOL0001.fasta','fasta'):
    seqs_list.append(record)

pro_list = []
for el in seqs_list:
    pro = el.translate()
    pro.id = el.id
    pro_list.append(pro)
    print(pro.id, len(pro.seq))
    sl(1)

pro_align = MultipleSeqAlignment(pro_list)

seq1 = SeqRecord(Seq('ATGTCTCGT', alphabet=generic_dna), id='pro1')
seq2 = SeqRecord(Seq('ATGCGT', alphabet=generic_dna), id='pro2')
seq3 = SeqRecord(Seq('ATGTCTCGT', alphabet=generic_dna), id='pro3')
s = [seq1,seq2,seq3]

pro1 = SeqRecord(Seq('MSR', alphabet=generic_protein), id='pro1')
pro2 = SeqRecord(Seq('M-R', alphabet=generic_protein), id='pro2')
pro3 = SeqRecord(Seq('MSR', alphabet=generic_protein), id='pro3')
p = [pro1,pro2,pro3]

aln = MultipleSeqAlignment(p)
codon_aln = build(aln, s)
print(codon_aln)
# print(pro2)
# x = Bio.codonalign.codonseq.cal_dn_ds(codon_aln)
# print(x)



if __name__ == '__main__':
    main()
