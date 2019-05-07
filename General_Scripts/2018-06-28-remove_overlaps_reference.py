from time import sleep as sl
import pandas as pd

infile = open('/Users/liam/Desktop/reference/GCA_000006745.gff', 'r').read().splitlines()

info_list = []
for line in infile:
    if 'CDS' in line and 'ID=' in line:
        col = line.split('\t')
        contig = col[0]
        gene_id = col[8].split(';')[0].split('_')[0][3:]

        beg_cord = col[3]
        end_cord = col[4]
        size = int(end_cord) - int(beg_cord)
        info_list.append(contig + '\t' + gene_id + '\t' + beg_cord + '\t' + end_cord + '\t' + str(size))

non_overlap = 0
nonlist = []
overlap = 0
overlist = []
equal = 0
equallist = []
cds_count = 0

outfile = open('/Users/liam/Desktop/reference_info.csv', 'w')
outfile.write('Gene' + ',' + 'Overlap' + ',' + 'Not Overlaping' + ',' + 'Equal' + ',' + 'Total Genes' + '\n')
for stat in info_list:
    col = stat.split('\t')
    #the second gene start
    if 'cds' + str(cds_count) == col[1]:
        gene = col[1]
        end = int(col[3])
        #gets the first gene end
        for nstat in info_list:
            ncol = nstat.split('\t')
            if 'cds' + str(cds_count + 1) == ncol[1]:
                start = int(ncol[2])
                #compares start and end of gene 1 and gene 2
                if start > end:
                    non_overlap = non_overlap + 1
                    nonlist.append(stat + '\t' + str(start) + '\t' + str(end))
                if start < end:
                    overlap = overlap + 1
                    overlist.append(stat + '\t' + str(start) + '\t' + str(end))
                if start == end:
                    equal = equal + 1
                    equallist.append(stat + '\t' + str(start) + '\t' + str(end))
        cds_count = cds_count + 1

for line in equallist:
    print(line)

# print('Overlap = ' + str(overlap))
# print('NonOverlap = ' + str(non_overlap))
# print('Equal = ' + str(equal))
# print('Total = ' + str(overlap + non_overlap + equal))
# print('Real Total = ' + str(len(info_list)))
