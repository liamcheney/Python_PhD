from time import sleep

infile =  open('/Users/liam/Desktop/testvis/clustered_proteins', 'r').readlines()
outfile = open('/Users/liam/Desktop/clustered_proteins1', 'w')

tag_list = []
temp_list = []
for i in infile:
    col = i.split('	')
    for j in col:
        col1 = j.split('::')[0]
        tag_list.append(str(col1))
    tag_list.append('\n')
    final = ' '.join(tag_list)

print(final)

outfile.write(final)
outfile.close()