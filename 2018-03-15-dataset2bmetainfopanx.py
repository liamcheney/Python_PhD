import json
import glob

strain_dict = {}
date_dict ={}
name_dict ={}
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
                date_dict[faccfile] = date
            elif '# Organism name:  ' in line:
                orgname = line.lstrip('# Organism name:  ')
                orgname1 = orgname.rstrip('(g-proteobacteria)\n')
                name_dict[faccfile] = orgname1

metafile = open('/Users/liam/pan-genome-visualization/public/dataset/dataset3eshort/strainMetainfo.json')
metain = json.load(metafile)
outfile = open('/Users/liam/pan-genome-visualization/public/dataset/dataset3eshort/realstrainMetainfo.json', 'w')

for key in metain:
    maininfo = metain[key]
    for i in maininfo:
        acc = i['accession']
        acc = acc.split('_031')[0]
        if "GCA_000006745" in acc:
            acc = acc.split(".")[0]
        i['collection_date'] = date_dict[acc]
        i['strain'] = strain_dict[acc]
        i['organism'] = name_dict[acc]
json.dump(maininfo, outfile)
