from Bio import SeqIO
import re

inseqs = SeqIO.parse("/Users/liamcheneyy/Desktop/Vibrio_D3A_alleles_ref.fasta","fasta")

for i in inseqs:
   s = str(i.seq)
   gene = str(i.id).split(":")[0]
   test = re.findall("(A{8,}|T{8,}|G{8,}|C{8,})",s)
   if len(test) >0:
       print(gene,"\t".join(test))
       # print(gene)