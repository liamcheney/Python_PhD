import glob

for i in range(14,15,1):
    print(i)
    gene_dict = {}
    for filename in glob.iglob('/Users/liamcheneyy/Desktop/vcseventh_' + str(i) + '/Mgtfi_ref/Schemes/*accessions*txt'):
        scheme_id = filename.split('/')[-1].split('_')[0]
        gene_dict[scheme_id] = {}
        infile = open(filename,'r').read().splitlines()
        for line in infile:
            gene_dict[scheme_id][line.strip()] = '1'

    print(gene_dict.keys())
    for key in gene_dict.keys():
        with open('/Users/liamcheneyy/Desktop/vcseventh_' + str(i) + '/Mgtfi_ref/refonly_allelic_profiles/' + key + '_gene_profiles.txt','w') as out:
            out.write('ST' + '\t' + 'dST' + '\t')
            for el in gene_dict[key]:
                out.write(el + '\t')
            out.write('\n')
            out.write('1' + '\t' + '0' + '\t')
            for val1 in gene_dict[key]:
                out.write(gene_dict[key][val1] + '\t')
