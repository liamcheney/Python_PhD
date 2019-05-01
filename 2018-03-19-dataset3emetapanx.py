import json
import glob

strain_dict = {}
date_dict = {}
name_dict = {}
accession_list = []

for filename in glob.glob('/Users/liam/Downloads/genome_assemblies/ncbi-genomes-2018-03-14/GCA*.txt'):
    accfile = filename.rstrip('_assembly_report.txt \n').replace('/Users/liam/Downloads/genome_assemblies/ncbi-genomes-2018-03-14/', '')
    faccfile = accfile[0:13]
    accession_list.append(faccfile)
    with open(filename, 'r') as infile:
        for line in infile:
            if '# Infraspecific name:  ' in line:
                strain = line.strip('# Infraspecific name:  strain= \n')
                strain_dict[faccfile] = strain
            elif '# Date:           ' in line:
                date = line.strip('# Date:           \n')
                date1 = date[0:4]
                date_dict[faccfile] = date1
            elif '# Organism name:  ' in line:
                orgname = line.lstrip('# Organism name:  ')
                orgname1 = orgname.rstrip('(g-proteobacteria)\n')
                name_dict[faccfile] = orgname1

infile1 = open('/Users/liam/pan-genome-visualization/public/dataset/dataset3eshort/coreGenomeTree.json','r').read()
infile2 = open('/Users/liam/Google Drive/Honours/3Vibrio/StrainsList.txt', 'r').readlines()
outfile1= open('/Users/liam/pan-genome-visualization/public/dataset/dataset3eshort/coreGenomeTree.json11', 'w')


node = infile1.split("name\": \"")
strlist = []
tempfile = []
for i in infile2:
    str = i[0:13]
    strlist.append(str)

for j in node:
    rowlab = j[0:13]
    if rowlab in strlist:
        nnode = j.replace('"strain": "strain"', '"strain": ' + '"' + strain_dict[rowlab] + '"').replace('{""collection_date": "unknown"', '"collection_date": "' + date_dict[rowlab] + '"')
        tempfile.append(nnode)
    else:
        oline = j
        tempfile.append(j)

final = 'name\": \"'.join(tempfile)
outfile1.write(final)
outfile1.close()




