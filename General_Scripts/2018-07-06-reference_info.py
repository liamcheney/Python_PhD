from time import sleep as sl

cds_in = list(open('/Users/liam/Desktop/gene_list.txt', 'r').read().splitlines())
igr_in = list(open('/Users/liam/Desktop/IGR_regions.txt','r').read().splitlines())
reference_in = open('/Users/liam/Desktop/reference/GCA_000006745.gff', 'r').read().splitlines()
outfile = open('/Users/liam/Desktop/D3A_gene_reference_info.csv', 'w')


def reference_information_dict_creator(reference_gff):
    info_dict = {}
    chromosome_list = []
    cds_lines = []
    for i in reference_gff:
        if 'CDS' in i and 'ID=' in i:
            cds_lines.append(i)
            col = i.split('\t')
            chromosome_list.append(col[0])
    chromosome_list = list(set(chromosome_list))

    gene_dict = {}
    cds_dict = {}
    for i in reference_gff:
        if 'ID=gene' in i:
            locus_tag = i.split('\t')[8].split(';')[-1].split('locus_tag=')[-1]
            gene = i.split('\t')[8].split(';')[0].split('_')[0].split('ID=')[-1]
            number = gene.split('gene')[-1]
            gene_dict[gene] = locus_tag
        if 'ID=cds' in i:
            cds = i.split('\t')[8].split(';')[0].split('_')[0].split('ID=')[-1]
            parent = i.split('\t')[8].split(';')[1].split('Parent=')[1]
            cds_dict[cds] = parent

    for i in range(len(cds_lines) - 1):
        line = cds_lines[i]
        col = line.split('\t')
        chromosome = col[0]
        start = col[3]
        end = col[4]
        strand = col[6]
        cds = col[8].split(';')[0].split('_')[0].split('ID=')[1]
        cds_num = int(cds.split('cds')[1])
        parent = col[8].split(';')[1].split('Parent=')[1]

        nextline = cds_lines[i + 1]
        nextcol = nextline.split('\t')
        nextstart = nextcol[3]
        nextend = nextcol[4]
        nextstrand = nextcol[6]
        nextcds = nextcol[8].split(';')[0].split('_')[0].split('ID=')[1]
        nextcds_num = int(nextcds.split('cds')[1])

        for chr in chromosome_list:
            if chr == chromosome:
                if end >= nextstart:
                    difference = int(end) - int(nextstart) + 1
                    newend = int(end) - difference
                    info_dict[cds] = {'start': start, 'end': newend, 'strand': strand, 'chromosome': chromosome, 'locustag': gene_dict[cds_dict[cds]]}

                    if difference == 0:
                        newend = int(end) - 1
                        info_dict[cds] = {'start': start, 'end': newend, 'strand': strand, 'chromosome': chromosome, 'locustag': gene_dict[cds_dict[cds]]}

                elif end >= nextstart and end >= nextend:
                    pass

                else:
                    info_dict[cds] = {'start': start, 'end': end, 'strand': strand, 'chromosome': chromosome, 'locustag': gene_dict[cds_dict[cds]]}
    return info_dict

def cds_list(reference_gff, cds_in, outfile):
    cds_list = []
    info_dict = reference_information_dict_creator(reference_gff)
    for key in info_dict.keys():
        if key in cds_in:
            outfile.write(info_dict[key]['locustag'] + ',' + str(info_dict[key]['start']) + ',' + str(info_dict[key]['end']) + ',' + info_dict[key]['strand'] + ',' + info_dict[key]['chromosome'] + '\n')
            outfile.write('\n')

    return cds_list

def igr_list(reference_gff, igr_in, outfile):
    igr_list = []
    info_dict = reference_information_dict_creator(reference_gff)
    for line in igr_in:
        col = line.split('_')
        firstcds = col[0]
        secondcds = col[1]
        outfile.write(info_dict[firstcds]['locustag'] + '_' + info_dict[secondcds]['locustag'] + ',' + str(int(info_dict[firstcds]['end']) + 1) + ',' + str(int(info_dict[secondcds]['start']) -1) + ',' + '+' + ',' + info_dict[firstcds]['chromosome'] + '\n')
    return igr_list

cds_list(reference_in, cds_in, outfile)
igr_list(reference_in, igr_in, outfile)