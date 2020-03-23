from time import sleep as sl

infile = open('/Users/liam/Desktop/testvis/clustered_proteins', 'r').readlines()
outfile = open('/Users/liam/d3e100clustered_proteins1', 'w')

group_list = []
tempfile =[]

for i in infile:
    line = i.split('\t')
    linelist = []
    for k in line:
        col = k.split('::')
        for l in col:
            if "GCA" in l:
                linelist.append(l)
            elif 'group' in l:
                linelist.append(l)
            elif 'cds' in l:
                linelist.append(l)

    final  = '\t'.join(linelist)
    tempfile.append(final)

final1 = '\n'.join(tempfile)
print(final1)
sl(1)
# outfile.write(str(final1))
# outfile.close()

# tempfile1 = []
#
# for i in infile:
#     group = i.split('group')
#     for j in group:
#         newline = j.split('\n')
#         for k in newline:
#             col = k.split('::')
#             print(col)
#             if "GCA" in col:
#                 tempfile1.append(col)
#                 print(col)
#             elif 'group' in col:
#                 tempfile1.append(col)
#
# print(tempfile1)