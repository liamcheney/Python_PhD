import glob
outfile = open('/Users/liam/Desktop/d3emetadata.tsv', 'w')

strainlist= []
acclist = []
datelist =[]
orgnamelist =[]

for filename in glob.glob('/Users/liam/Downloads/genome_assemblies/ncbi-genomes-2018-03-14/GCA*.txt'):
    with open(filename, 'r') as infile:
        for line in infile:
            if '# GenBank assembly accession: ' in line:
                acc = line.replace('# GenBank assembly accession: ', '')[0:13]
                acclist.append(acc)
            elif '# Infraspecific name:  ' in line:
                strain = line.replace('# Infraspecific name:  strain=','').rstrip("\n")
                strainlist.append(strain)
            elif 'Isolate' in line:
                print(line)
            elif '# Date:           ' in line:
                date = line.strip('# Date:           \n')
                datelist.append(date)
            elif '# Organism name:  ' in line:
                orgname = line.replace('# Organism name:  ','').rstrip('(g-proteobacteria)\n')
                orgnamelist.append(orgname)
