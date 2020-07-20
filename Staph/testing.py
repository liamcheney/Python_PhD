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
# save_list = []
# for records in SeqIO.parse('/Users/liamcheney/Desktop/refs.fasta','fasta'):
#     locus=records.id.split(':')[0]
#     if locus not in save_list:
#         save_list.append(locus)
#
# for element in save_list:
#     lists = []
#     for record in SeqIO.parse('/Users/liamcheney/Desktop/refs.fasta', 'fasta'):
#         in_locus=record.id.split(':')[0]
#         if element == in_locus:
#             lists.append(record)
#     SeqIO.write(lists, '/Users/liamcheney/Desktop/locus/' + element + '.fasta','fasta')
#     print(element, len(lists))

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


    infile = open('/Users/liamcheney/Desktop/meta_7-gene-alleles.txt').read().splitlines()
    save = {}
    for line in infile[1:]:
        col = line.split('\t')
        biopro = col[3]
        if biopro not in save.keys():
            save[biopro] = 0
    fin = {}
    listss = []
    for i in save.keys():
        count  = 0
        both_meta = 0
        for line in infile[1:]:
            col = line.split('\t')
            biopro = col[3]
            year = col[4]
            country = col[5]

            if i == biopro:
                count += 1

                if year != '' and country != '':
                    both_meta += 1



        fin[i] = count
        print(i, count, both_meta, int(both_meta/count * 100))

    with open('/Users/liamcheney/Desktop/test-gene-alleles.txt','w') as out:
        for el in listss:
            for i in el:
                out.write(i + '\t')
            out.write('\n')
            out.write(el[0])


if __name__ == '__main__':
    main()
