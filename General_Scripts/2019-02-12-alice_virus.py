from time import sleep as sl
import glob
from Bio import SeqIO
from Bio.Seq import Seq

input = open("/Users/liamcheneyy/Desktop/isca_NIRVS_positions.csv", 'r').read().splitlines()

info_dict = {}
count = 0
for line in input:
    count = count + 1
    col = line.split(',')
    scaff_num = col[0].strip('\ufeff').strip() + "." + str(count)
    start = col[1]
    end = col[2]
    info_dict[scaff_num] = {"start":start, "end": end}

fasta_in = "/Users/liamcheneyy/Desktop/Ixodes_scapularis.IscaW1.dna.toplevel.fasta"
fast_dict = {}
with open(fasta_in) as file:
    for sequence in SeqIO.parse(file, 'fasta'):
        scaff_num_1 = sequence.id
        seq_1 = sequence.seq
        fast_dict[scaff_num_1] = str(seq_1)

out = open("/Users/liamcheneyy/Desktop/output_Ixodes_scapularis.IscaW1.dna.toplevel.fasta", 'w')
for key, value in info_dict.items():
    seq = fast_dict[key.split('.')[0]]
    start = int(value["start"]) - 1 - 5000
    end = int(value["end"]) - 1 + 5000
    want_seq = seq[start:end]
    out.write(">" + key + '\n')
    out.write(want_seq + '\n')

out.write(">DS971868.142" + '\n')
out.write(fast_dict["DS971868"] + '\n')

out.close()