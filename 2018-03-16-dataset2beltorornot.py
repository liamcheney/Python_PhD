import json
infile2 = open('/Users/liam/Desktop/7pangenlist1.txt', 'r').readlines()
infile1 = open('/Users/liam/pan-genome-visualization/public/dataset/dataset2bshort/coreGenomeTree.json','r').read()
outfile = open('/Users/liam/pan-genome-visualization/public/dataset/dataset2bshort/coreGenomeTree1.json', 'w')

node = infile1.split('{"name": "GCF_')

list = []
list2 = []

for j in infile2:
    var = j.strip('\n')
    list.append(var)

#will go through each list inside node. spliting creates lists. each list is added in order
#based on presence of collection date or not.
for i in node:
    rowlab = i[0:9]
    if rowlab in list:
        newl = i.replace('{"collection_date": "unknown"', '{"collection_date": "7th pandemic"')
        list2.append(newl)
    else:
        oline = i
        list2.append(oline)

#adds the identifier i split by before and creates the file as normal before.
nfile = '{"name": "GCF_'.join(list2)
outfile.write(nfile)
outfile.close()



