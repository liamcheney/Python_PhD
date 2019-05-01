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
                strain = line.replace('# Infraspecific name:  strain=', '').replace('\n','')
                strain_dict[faccfile] = strain
            elif '# Date:           ' in line:
                date = line.replace('# Date:           ','').replace('\n', '')[0:4]
                date_dict[faccfile] = date
            elif '# Organism name:  ' in line:
                orgname = line.lstrip('# Organism name:  ')
                orgname1 = orgname.rstrip('(g-proteobacteria)\n')
                name_dict[faccfile] = orgname1

infile1 = open('/Volumes/Macintosh HD/Users/liam/d3eStrainsList.txt', 'r').readlines()
panlist = []

for i in infile1:
    panacc = i[0:13]
    panlist.append(panacc)

for j in panlist:
        print(name_dict[j])
#
# for k in panlist:
#     print(date_dict[k])
# #
# for l in panlist:
#     print(strain_dict[l])