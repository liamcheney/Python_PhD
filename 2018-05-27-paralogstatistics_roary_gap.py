from time import sleep as sl
import math
import pandas
from time import clock

print(str(clock()) + '\t' + 'Beginning')
#overall the script determines the number of paralog genes for each strains in roary analysis and also relates them to genome lengths.
infile = open('/Users/liam/Desktop/gene_presence_absence.csv','r').readlines()
outfile1 = open('/Users/liam/Desktop/kept_lines.csv','w')

strains_used = 0
print(str(clock()) + '\t' + 'caculating paralogs groups in each ortholog')
outfile1.write(infile[0])
#creates a dictionary with the number of cells containing genes for each ortholog group
for i in infile[2:]:
    cell_count_with = 0
    cell_count_empty = 0
    col = i.split(',')
    strains_used = len(col[14:])
    for j in col[14:]:
        if '_' in j:
            cell_count_with += 1
        else:
            cell_count_empty += 1

        total = cell_count_empty + cell_count_with
        if total == strains_used:
            if cell_count_with >= math.ceil(strains_used * 0.99):
                outfile1.write(str(i))

print(str(clock()) + '\t' + 'creating dictionary with paralogs in each ortholog')
#iterates over each column (strain) and counts the number of orthologs which are paralogous. Dictionary is created with strain:paralogous group count.
kept = open('/Users/liam/Desktop/kept_lines.csv', 'r').read().splitlines()

strain_orth_dict ={}
strains_names = kept[0].split(',')
for t in strains_names[14:]:
    strains_strip = t.lstrip('"').rstrip('"\n')
    strain_orth_dict[strains_strip] = 'x'

new_dict_keys = {}
for j in strain_orth_dict.keys():
    if "GCF" in j:
        new_key = "GCA" + j[3:].rstrip('"\n')
        new_dict_keys[new_key] = 'x'
    else:
        new_dict_keys[j] = 'x'

print(str(clock()) + '\t' + 'calculating number of paralogs in each strain')
#counts number of paralog groups in each strain and the strain name for each list. each list is then combined into final dict.
paralog_count_list = [0 for x in range(strains_used)]
paralog_group_list = [0 for x in range(strains_used)]
for j in kept[1:]:
    col2 = j.split(',')[14:]
    orth_name = (col2[0].lstrip('"').rstrip('"'))
    for c in range(strains_used):
        paralog_count = str(col2[c]).count('\t') + 1
        paralog_count_list[c] = paralog_count
        if paralog_count_list[c] > 1:
            paralog_group_list[c] = paralog_group_list[c] + 1

final_dict=dict(zip(new_dict_keys.keys(), paralog_group_list))

# Only needed the first time get stats, were all written out to lengths file.
print(str(clock()) + '\t' + 'extracting NCBI and raw data assembly statistics')
#extracts the genome statistics for all the NCBI assembled strains
import glob

NCBI_length_dict = {}
for filepath in glob.glob('/Users/liam/Desktop/genome_assemblies/ncbi-genomes-2018-05-24/*'):
    with open(filepath,'r') as file:
        filename = filepath.split('/')[6][0:13]
        for i in file:
            if 'all	all	all	all	total-length	' in i:
                NCBI_length_dict[filename] = int(i.split('\t')[5].rstrip('\n'))

#creates dictionary of keys = strains accession and values = their total genome length
infile1 = open('/Users/liam/Desktop/lengths.csv','r').readlines()
genome_size_dict = {}
for t in infile1[1:]:
    col1 = t.split(',')
    genome_size_dict[col1[0]] = col1[1].lstrip('"').rstrip('\n')

#reads in a list of 7th pandemic strains. used for isolating 7th pandemic from D2C pandemic.
seventh_in = open('/Users/liam/Desktop/seventh-strains.txt', 'r').readlines()
seventh_list = []
for i in seventh_in:
    seventh_list.append(i.rstrip('\n'))

#reads in list of non-seventh pandemic from the D2C.
not_seventh = open('/Users/liam/Desktop/not seventh.txt', 'r').readlines()
not_seven = []
for i in not_seventh:
    not_seven.append(i.rstrip('\n'))

print(str(clock()) + '\t' + 'writing out orthologo groups with no paralogs')
#print final list of paralogs groups and length of each genome
outfile2 = open('/Users/liam/Desktop/D3A_genome_length_vs_paralogs.csv','w')
outfile2.write('Pandemic' + ',' + 'Accession' + ',' + 'Length' + ',' + 'Paralog Groups' + '\n')

GCA_list = final_dict.keys()
count_sev = 0
not_count = 0
for u in GCA_list:
    if u in seventh_list:
        count_sev = count_sev + 1
        outfile2.write('Seventh' + ',' +  str(u) + ',' + str(genome_size_dict.get(u, 'None')) + ',' + str(final_dict.get(u, 'None')) + '\n')
    elif u in not_seven:
        not_count = not_count + 1
        outfile2.write('Not Sevnth' + ',' + str(u) + ',' + str(genome_size_dict.get(u, 'None')) + ',' + str(final_dict.get(u, 'None')) + '\n')
outfile2.close()

print(count_sev)
print(not_count)