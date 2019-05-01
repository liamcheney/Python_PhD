from time import sleep as sl
from Bio import SeqIO
import glob

ref_alleles = SeqIO.parse('/Users/liam/Desktop/Vibrio_D3A_alleles_ref.fasta', 'fasta')
reference_in = SeqIO.parse('/Users/liam/Desktop/Reference Info/GCF_000006745.1_ASM674v1_genomic.fna', 'fasta')
out_file = open('/Users/liam/Desktop/Vibrio_D3A_alleles_ref_test.fasta','w')

chrI = ''
chrII = ''
for sequence in reference_in:
    if sequence.id == 'NC_002505.1':
        chrI = str(sequence.seq)
    if sequence.id == 'NC_002506.1':
        chrII = str(sequence.seq)

for allele in ref_alleles:
    sequence1 = str(allele.seq)
    if (sequence1 in chrI) and (sequence1 in chrII):
        out_file.write(allele.id[:-2] + '\t' + "Both" + '\n')
    elif sequence1 in chrI:
        out_file.write(allele.id[:-2] + '\t' + "ChrI" + '\n')
    elif sequence1 in chrII:
        out_file.write(allele.id[:-2] + '\t' + "ChrII" + '\n')
    elif (sequence1 not in chrI) and (sequence1 not in chrII):
        out_file.write(allele.id[:-2] + '\t' + "None Found" + '\n')

out_file.close()