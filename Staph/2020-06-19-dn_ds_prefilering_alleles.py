import Bio.SeqIO
import glob
import numpy as np

infile = open('/Users/liamcheney/Desktop/Book2.txt', 'r').read().splitlines()
size_dict = {}
for line in infile:
    col = line.split('\t')
    size_dict[col[0]] = int(col[1])

save_dict = {}
for filename in glob.iglob('/Users/liamcheney/Desktop/filtered_alleles/SACOL*.fasta'):
    save_list = []
    locus = filename.split('/')[-1].strip('.fasta')

    for records in Bio.SeqIO.parse(filename, 'fasta'):
        if '-' not in records.id:
            length = len(records.seq)
            save_list.append(length)
            # if length == size_dict[locus]:
            #     save_list.append(records)

    save_dict[locus] = save_list

    # for k,v in save_dict.items():
    #     Bio.SeqIO.write(v, '/Users/liamcheney/Desktop/filtered_alleles/' + k + '.fasta','fasta')

    freq = {x: save_list.count(x) for x in save_list}
    save_dict[locus] = freq

for k, v in save_dict.items():
    total = sum(v.values())
    largest = list(v.keys())[0]
    num_w_larg = v[list(v.keys())[0]]
    coveraged = round(num_w_larg / total * 100, 2)

    print(k, len(v), total, largest, coveraged, sep='\t')