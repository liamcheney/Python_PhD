import glob
import os
import shutil

# in_assem = open('/Users/liam/Desktop/D2C_list.txt', 'r').read().splitlines()
#
# used_assemblies = [i for i in in_assem]
#
# for i in glob.iglob(r'/Users/liam/Desktop/ncbi-genomes-2018-06-07/*'):
#     with open(i, 'r') as file:
#         filename = i.split('/')[-1][0:13]
#         for j in used_assemblies:
#             if j in filename:
#                 print(j)
#                 shutil.move(i, "/Users/liam/Desktop/keptstrains/")

NCBI_gene_count_dict = {}
for i in glob.iglob('/Users/liam/Desktop/keptstrains/*'):
    with open(i, 'r') as file:
        gene_count = 0
        filename = i.split('/')[-1][0:13].replace('GCF','GCA')
        for line in file:
            if '#' not in line and'protein_coding' in line:
                col = line.split('\t')
                if 'gene' == col[2]:
                    gene_count = gene_count + 1
    NCBI_gene_count_dict[filename] = gene_count

prokka_gene_count_dict = {}
for i in glob.iglob('/Users/liam/Desktop/D2C_prokka_gffs/*'):
    filename = i.split('/')[-1][0:13]
    if filename in list(NCBI_gene_count_dict.keys()):
        with open(i, 'r') as file:
            infile = file.read().splitlines()
            gene_count = 0
            filename = i.split('/')[-1][0:13]
            for line in infile:
                if 'Prodigal:2.6' in line:
                    col = line.split('\t')
                    if 'CDS' == col[2]:
                        gene_count = gene_count + 1
        prokka_gene_count_dict[filename] = gene_count

print(len(NCBI_gene_count_dict.keys()))
print(len(prokka_gene_count_dict.keys()))

outfile = open('/Users/liam/Desktop/NCBI_vs_Prokka_gene_counts.csv', 'w')

outfile.write('Accession' + ',' + 'Gene Count' + '\n')
for i,x in NCBI_gene_count_dict.items():
    outfile.write('NCBI' + ',' + str(i) + ',' + str(x) + '\n')

for k,n in prokka_gene_count_dict.items():
    outfile.write('Prokka' + ',' + str(k) + ',' + str(n) + '\n')

