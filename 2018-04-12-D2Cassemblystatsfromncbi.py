import glob2
import numpy

outfile = open('/Users/liam/Desktop/D2Cassembly_coverages.csv', 'w')

nostats = []
for filename in glob2.iglob('/Users/liam/Desktop/genome_assemblies/ncbi-genomes-2018-04-11/*stats*'):
    acc = filename.replace('/Users/liam/Desktop/genome_assemblies/ncbi-genomes-2018-04-11/','')[0:13]
    infile = open(filename, 'r').read()
    # print(acc)
    # for i in infile:
    #     if '# Genome coverage:' in i:
    #         col = i.split('# Genome coverage: ')
    #         cov = str(col[1]).rstrip('x\n').rstrip('X\n')
    #         print(str(acc) + '\t' + str(cov))
    if '# Genome coverage: ' not in infile:
        nostats.append(filename)

for i in nostats:
    infile2 = open(i, 'r').readlines()
    acc = i.replace('/Users/liam/Desktop/genome_assemblies/ncbi-genomes-2018-04-11/','')[0:13]
    for j in infile2:
        if '# Assembly level: Complete ' in j:
            print(acc)
