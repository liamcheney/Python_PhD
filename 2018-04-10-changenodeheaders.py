import sys
from time import sleep as sl

input = sys.argv[1]
filename = input[0:-4]

infile = open(input,'r').readlines()
outfile = open(input,'w')

for i in infile:
    if '>NODE' in i:
        col = i.split('_')
        col[0] = '>' + filename + '_NODE'
        new = col[0:2]
        newr = '_'.join(new)
        outfile.write(newr + '\n')
    else:
        outfile.write(i)

outfile.close()