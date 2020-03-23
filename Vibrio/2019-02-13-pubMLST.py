# from time import sleep as sl
# infile = open('/Users/liamcheneyy/Desktop/results.tsv', 'r').read().splitlines()
# table = open('/Users/liamcheneyy/Desktop/GCF_000006745.1_ASM674v1_feature_table.txt').read().splitlines()
#
# for line in infile:
#     vc = ""
#     col = line.split('\t')
#     NP = col[0].split('_')[-3] + '_' + col[0].split('_')[-2]
#     for line1 in table:
#         if NP in line1:
#             col1 = line1.split('\t')
#             vc = col1[-4]
#     print(line + '\t' + vc + '\t' + NP)

# gff = open('/Users/liamcheneyy/Desktop/GCA_000006745.gff').read().splitlines()
# core_genes = open('/Users/liamcheneyy/Desktop/D2C_core_genes.txt').read().splitlines()
#
# for line in core_genes:
#     cds = "cds" + line.split('_')[2]
#     for line1 in gff:
#         if "CDS" in line1:
#             col = line1.split(';')[0]
#             gff_cds = col.split('\t')[-1].split('_')[0][3:]
#             if cds == gff_cds:
#                 protein = line1.split(';')[2].split(':')[1].split(',')[0]
#                 print(protein)

from Bio import SeqIO
from Bio.Seq import Seq
from time import sleep as sl

##making fasta of nucl for unique MGT genes
# gff = '/Users/liamcheneyy/Desktop/GCF_000006745.1_ASM674v1_cds_from_genomic.fna'
# core_genes = open('/Users/liamcheneyy/Desktop/MGT_unique.txt').read().splitlines()
# outfile = open('/Users/liamcheneyy/Desktop/MGT_unique_seqs.txt','w')
#
# for line in core_genes:
#     print(line)
#     with open(gff,'r') as file:
#         for record in SeqIO.parse(file, 'fasta'):
#             if line in record.description:
#                 outfile.write('>' + record.description + '\n')
#                 outfile.write(str(record.seq) + '\n')
# outfile.close()

##making fasta of nucl for unique pubMLST genes
# blast_in = open('/Users/liamcheneyy/Desktop/pubmlst_blast.tsv', 'r').read().splitlines()
# core_genes = open('/Users/liamcheneyy/Desktop/pubmslst_unique.txt').read().splitlines()
# outfile = open('/Users/liamcheneyy/Desktop/pubMLST_unique_seqs.txt','w')
#
# core_genes_list = [x for x in core_genes]
# for line in blast_in:
#     col = line.split('\t')
#     vc = col[12]
#     vchol = col[1].split(':')[-1]
#     if vc in core_genes_list:
#         print(vc)
#         with open('/Users/liamcheneyy/Desktop/vc_pubmlst_core_genes.xmfa','r') as file:
#             for record in SeqIO.parse(file,'fasta'):
#                 if vchol in record.id:
#                     outfile.write('>' + record.description + '\n')
#                     outfile.write(str(record.seq) + '\n')

# converting protein names to cds names
# gff = open('/Users/liamcheneyy/Desktop/GCA_000006745.gff', 'r').read().splitlines()
# genes = open('/Users/liamcheneyy/Desktop/only_pub.txt','r').read().splitlines()
# genes_list = [x for x in genes]
#
# for line in gff:
#     if "CDS" in line:
#         col = line.split(';')
#         pro = col[2].split(',')[0].split(':')[-1]
#         cds = col[0].split('ID=')[-1].split('_')[0]
#         # if pro in genes_list:
#         #     print(pro + '\t' + cds)
#         if pro not in genes_list:
#             print(pro)

# #finding if any of the unique genes were found split
# roary_in = open('/Users/liamcheneyy/Desktop/96_gap.csv', 'r').read().splitlines()
# genes = open('/Users/liamcheneyy/Desktop/only_pubmlst_cds','r').read().splitlines()
# outfile = open('/Users/liamcheneyy/Desktop/roary_gap_pubmlst_genes.csv' ,'w')
# gene_list= [x for x in genes]
#
# count = 0
# for line in roary_in:
#     presence_count = 0
#     col = line.split(',')
#     for j in col[15:]:
#         if 'cds' in j:
#             for element in gene_list:
#                 if element + '_' in j:
#                     # outfile.write(line + '\n')
#                     if count("cds") > 2:
#                         count = count + 1
#
# outfile.close()
# print(count)
# percent = count / len(gene_list)
# print(str(percent))

# infile = open('/Users/liamcheneyy/Desktop/roary_gap_pubmlst_genes.csv','r').read().splitlines()
#
# for line in infile:
#     cut_off = 'No'
#     par_count = 0
#     cell_count = 0
#     cds = ""
#     col = line.split(',')
#     for j in col[15:]:
#         if '\t' in j:
#             par_count = par_count + 1
#         if '_' in j:
#             cell_count = cell_count + 1
#         if "cds" in j:
#             if '\t' in j:
#                 cds = j.split('\t')
#                 cds = ' '.join(cds)
#             else:
#                 cds = j.strip('"')
#     if cell_count >= (0.99 * 2282):
#         cut_off = "Yes"
#     print(cds + '\t' + str(cell_count) + '\t' + str(par_count) + '\t' + cut_off)

