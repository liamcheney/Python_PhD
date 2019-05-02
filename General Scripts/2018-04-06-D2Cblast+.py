# import glob2
# panlist = []
# notlist = []
# for filename in glob2.iglob('/Users/liam/Desktop/Database/DB*/*.fasta'):
#     db = filename.replace('/Users/liam/Desktop/Database/','').rstrip('/.fasta')
#     infile = open(filename, 'r').readlines()
#     for line in infile:
#         if 'ERR' in line:
#             col = line.split('\t')
#             if col[3] == '100.000':
#                 panlist.append(line)
#             if col[3] != '100.000':
#                 notlist.append(line)
#
# for j in panlist:
#     col1 = j.split('\t')
#     err = col1[0]
#     err = err[0:9]
#     allele = col1[1]
#     percen = col1[3]
#     # print(err + '\t' + allele + '\t' + percen)
#     print(err)