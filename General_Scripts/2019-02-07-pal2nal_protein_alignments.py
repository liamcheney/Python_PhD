from time import sleep as sl
import glob
from Bio import SeqIO
from Bio.Seq import Seq
import progressbar
from time import sleep
import sys

# #checking stats of pal2nal oputputs, based on muscle alignments
# count_list = []
# genomes_used = 0
#
# for filename in glob.iglob('/srv/scratch/z5087966/HGT_seventh_pandemic/2019-01-30-core_genes_DNDS/Pal2nal/*fna'):
#     ATG_count = 0
#     GTG_count = 0
#     TTG_count = 0
#     out_name = filename.split('/')[-1][:-12]
#     print(out_name)
#     with open(filename) as file:
#         for record in SeqIO.parse(file, 'fasta'):
#             genomes_used = genomes_used + 1
#             if record.seq[0:3] == 'ATG':
#                 ATG_count = ATG_count + 1
#             if record.seq[0:3] == 'GTG':
#                 GTG_count = GTG_count + 1
#             if record.seq[0:3] == 'TTG':
#                 TTG_count = TTG_count + 1
#     count_list.append([out_name, ATG_count, GTG_count, TTG_count, genomes_used])
#
# outfile = open('/srv/scratch/z5087966/HGT_seventh_pandemic/2019-01-30-core_genes_DNDS/Pal2nal/pal2nal_ATG_stats.csv', 'w')
# outfile.write('Group Name' + ',' + 'ATG Starts' + ',' + 'GTG Starts' + ',' + 'TTG Starts' + ',' + 'Number of Sequences' + '\n')
# for element in count_list:
#     for j in element:
#         outfile.write(str(j) + ',')
#     outfile.write('\n')

# outfile.close()

##checking the quality of input nucleotide fasta and protein
# count1_list = []
# for filename in glob.iglob("/Volumes/Liam's HDD/PhD/2_HGT_seventh_pandemic/2019-01-30-core_genes_DNDS/MSA_fas/*fna"):
#     N_count = 0
#     out_name = filename.split('/')[-1][:-4]
#     print(out_name)
#     with open(filename) as file:
#         for record in SeqIO.parse(file, 'fasta'):
#             if 'N' in record.seq:
#                 N_count = N_count + 1
#     count1_list.append([out_name, N_count])
#
#
# outfile = open('/Users/liamcheneyy/Desktop/MSA_fas_stats.csv', 'w')
# outfile.write('Group Name' + ',' + 'Number of Seqs with N' + '\n')
# for element in count1_list:
#     for j in element:
#         outfile.write(str(j) + ',')
#     outfile.write('\n')
#
# outfile.close()

# # writing out protein alignments
# for filename in glob.iglob('/Users/liamcheneyy/Desktop/subset_50_genomes/3_MSA_fna/*fna'):
#     out_name = filename.split('/')[-1][:-4]
#     with open(filename) as file:
#         out = open('/Users/liamcheneyy/Desktop/subset_50_genomes/4_MSA_aln/' + out_name + '.aln', 'w')
#         for seq in SeqIO.parse(file, 'fasta'):
#             prot = Seq.translate(seq.seq)
#             out.write('>' + seq.id + '\n')
#             out.write(str(prot) + '\n')
#     print(out_name)
#             # out.write(seq.id[0:10] + '\t' + str(seq.seq) + '\n')

pal2nal_fna = sys.argv[1]
output_folder = sys.argv[2]

##writing out codon alignments in SNAP tsv format
for filename in glob.iglob(pal2nal_fna + '/*'):
    out_name = filename.split('/')[-1][:-12]
    with open(filename) as file:
        print(out_name)
        out = open(output_folder + '/' + out_name + '.fna', 'w')
        for seq in SeqIO.parse(file, 'fasta'):
            out.write(seq.id[0:10] + '\t' + str(seq.seq) + '\n')
        out.close()
