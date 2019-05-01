from time import sleep as sl

infiel1 = open('/Users/liam/Desktop/gene_presence_absence.csv', 'r').readlines()
infile2 = open('/Users/liam/Desktop/test_gene_presence_absence.csv', 'r').readlines()
outfile = open('/Users/liam/Desktop/test_gene_presence_absence.csv', 'w')
outfile2 = open('/Users/liam/Desktop/strainsof_gene_presence_absence.csv', 'w')
instrain = 'ERR018111'

for i in infiel1:
    if "group_" in i:
        col = i.split(',')
        orthg = int(col[0].replace('"group_','').replace('"',''))
        seqnum = int(col[3].replace('"','').replace('"',''))
        avseq = float(col[5].replace('"','').replace('"',''))
        #writes the strains within range with certain avgseq to file.
        if 1751 <= seqnum <= 1769 and avseq == 1:
            outfile.write(i)

for j in infile2:
    if 'ERR018111' in col[14]:
        print(col[14])

outfile.close()

# creates a list of the strains in the .csv file. doesnt fully work for writing out.
    if "Gene" in i:
        col1 = i.split(',')
        for j in col1:
            if "ERR" in j:
                var = j.replace('','')
                print(j[1:10])
                outfile2.write(str(j[0:10]) + '\n')
            elif "GCA" in j:
                print(j[1:14])
                outfile2.write(str(j[0:14]) + '\n')
            elif "GCF" in j:
                print(j[1:14])
                outfile2.write(str(j[0:14]) + '\n')
    outfile2.close()
