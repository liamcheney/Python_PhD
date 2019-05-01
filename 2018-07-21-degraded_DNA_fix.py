from Bio import SeqIO
from Bio.Seq import Seq
import sys
infile = sys.argv[1]

fasta_sequence = SeqIO.parse(open(infile), 'fasta')

alleles_list = []
for fasta in fasta_sequence:
    new = fasta
    new.seq = Seq(str(fasta.seq).replace('R', 'N').replace('Y', 'N').replace('S', 'N').replace('W', 'N').replace('K', 'N').replace('M', 'N').replace('B', 'N').replace('D', 'N').replace('H', 'N').replace('V', 'N'))
    alleles_list.append(new)
SeqIO.write(alleles_list, sys.argv[2], 'fasta')






