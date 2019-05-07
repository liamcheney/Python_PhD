import sys
from time import sleep as sl

infile = open('/Users/liamcheneyy/Desktop/untitled folder/ERR018132.gbf','r').readlines()
outfile = open('/Users/liamcheneyy/Desktop/untitled folder/tERR018132.gbf','w')

list1 = []
for i in infile:
    if 'LOCUS' in i:
        col = i.split()
        # print(col)
        list = []
        if len(col) <= 6:
            col1 = col[1].split('_')
            # print(col1)
            bpl = -1*len(col[3])
            bp = col1[3]
            p1 = col1[5][:bpl]
            p12 = str(p1 + '    ' + bp)
            col1[5] = p12
            reb = '_'.join(col1)
            col[1] = reb
            reb1 = ' '.join(col)
            print(reb1)
            # outfile.write(str(reb1))
    #     else:
    #         reb = ' '.join(col)
    #         outfile.write(str(reb + '\n'))
    # else:
    #     outfile.write(i)
