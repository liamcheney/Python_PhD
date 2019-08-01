from time import sleep as sl

# cds_in = list(open('/Users/liamcheneyy/Desktop/gene_list.txt', 'r').read().splitlines())
igr_in = list(open('/Users/liamcheneyy/Desktop/ig_gene_list.txt','r').read().splitlines())
reference_in = open('/Users/liamcheneyy/Desktop/ref/GCF_000006745.1_ASM674v1_feature_table.txt', 'r').read().splitlines()
outfile = open('/Users/liamcheneyy/Desktop/D3A_gene_reference_info.csv', 'w')


def reference_information_dict_creator(reference_in):
    info_dict = {}
    chromosome_list = []
    cds_lines = []

    for i in reference_in:
        col = i.split('\t')
        if 'gene' in col[0]:
            locus_tag = col[16]
            start = col[7]
            end = col[8]
            chr = col[5]
            strand = col[9]
            length = int(end) - int(start)

            info_dict[locus_tag] = {'start': start, 'end': end, 'strand': strand, 'chromosome': chr, 'length':length}

    return info_dict

def cds_list(reference_in, locus_in, outfile):
    cds_list = []
    info_dict = reference_information_dict_creator(reference_in)
    for key in info_dict.keys():
        if key in locus_in:
            print(key)
            outfile.write(key + ',' + str(info_dict[key]['start']) + ',' + str(info_dict[key]['end']) + ',' + info_dict[key]['strand'] + ',' + info_dict[key]['chromosome'] + ',' + str(info_dict[key]['length']))
            outfile.write('\n')

    return cds_list

def igr_list(reference_in, igr_in, outfile):
    igr_list = []
    info_dict = reference_information_dict_creator(reference_in)
    for line in igr_in:
        col = line.split('_')
        firstcds = col[0]
        secondcds = col[1]
        start = str(int(info_dict[firstcds]['end']))
        end = str(int(info_dict[secondcds]['start']))
        length = int(end) - int(start)

        outfile.write(line.strip('\n') + ',' + start + ',' + end + ',' + str(length) + '\n')
    return igr_list

# cds_list(reference_in, cds_in, outfile)
igr_list(reference_in, igr_in, outfile)