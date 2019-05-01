import glob2
import csv
import pandas

outfile = open('/Users/liam/Desktop/d2Ccombinedstats.csv', 'w')

contlist = []
lenlist = []
L75 = []
N75 = []
L50 = []
N50 =[]
gene = []
geneunique = []
combinedlist = ('Accession' + '\t' + 'Number of Contigs' + '\t' + 'Total length' + '\t' + 'N50' + '\t' + 'L50' + '\t' + 'L75' + '\t' + 'N75' + '\t' + 'Predicted Genes' + '\t' + 'Unique Genes' + '\n')
outfile.write(combinedlist)

for filename in glob2.iglob('/Users/liam/Desktop/untitled folder 2/*'):
    with open(filename, 'r') as infile:
        fl = filename.replace('/Users/liam/Desktop/untitled folder 2/', '').replace('-report.txt', '')
        for line in infile:
            if '# contigs                   ' in line:
                cont = line.replace('# contigs                   ', '').replace('\n','')
                contlist.append(int(cont))
            if 'Total length                ' in line:
                len = line.replace('Total length                ', '').replace('\n','')
                lenlist.append(len)
            if 'N50                         ' in line:
                n50 = line.replace('N50                         ', '').replace('\n','')
                N50.append(n50)
            if 'L50                         ' in line:
                l50 = line.replace('L50                         ', '').replace('\n','')
                L50.append(l50)
            if 'L75                         ' in line:
                l75 = line.replace('L75                         ', '').replace('\n','')
                L75.append(l75)
            if 'N75                         ' in line:
                n75 = line.replace('N75                         ', '').replace('\n','')
                N75.append(n75)
            if '# predicted genes (>= 0 bp)' in line:
                gen = line.replace('# predicted genes (>= 0 bp)', '').replace('\n','')
                gene.append(gen)
            if '# predicted genes (unique)      ' in line:
                genuni = line.replace('# predicted genes (unique)      ', '').replace('\n','')
                geneunique.append(genuni)

    outfile.write(fl + '\t' + cont + '\t' + len + '\t' + n50 + '\t' + l50 + '\t' + n75 + '\t' + l75 + '\t' + gen + '\t' + genuni + '\t' + '\n')
outfile.close()



