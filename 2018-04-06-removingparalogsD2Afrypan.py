from time import sleep as sl

# infile1 = open('/Users/liam/FriPan/test2.proteinortho','r').readlines()
# outfile1 = open('/Users/liam/Desktop/test2.proteinortho','w')
#
# for i in infile1:
#     col = i.split('\t')
#     linelist = []
#     for j in col:
#         if ',' in j:
#             par = j.split(',')[0]
#             linelist.append(par)
#         else:
#             linelist.append(j)
#     final = '\t'.join(linelist)
#
#     outfile1.write(final)
#
# outfile1.close()

infile1 = open('/Users/liam/FriPan/test2.descriptions', 'r').readlines()
outfile1 = open('/Users/liam/Desktop/test2.descriptions', 'w')

for i in infile1:
    col = i.split(',')[0].strip('\n').strip('/t')
    outfile1.write(col + '\n')
outfile1.close()