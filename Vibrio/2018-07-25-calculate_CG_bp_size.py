from Bio.Seq import Seq
import sys

# infile_list= list(open(sys.argv[1]).read().splitlines())
# reference_gff_info = open(sys.argv[2], 'r').read().splitlines()
# reference_genome_length = sys.argv[3]

infile_list= list(open('/Users/liamcheneyy/Desktop/MGT9_gene_accessions.txt').read().splitlines())
reference_gff_info = open('/Users/liamcheneyy/Google Drive/Honours/Reference Info/gene_cds_locus.csv', 'r').read().splitlines()
reference_genome_length = 4033464


def calc_CG_bp_size(genes_list, reference_info, reference_genome_length):
    """
    :param genes_list: list of genes in the core genome, which will be summed to find total length.
    :param reference_info: a tsv file with the format of gene_id : {start:XX, end:XX, strand:XX, chromsome:XX}
    :param reference_genome_length: the total size of the reference genome (eg. 4100000 million b.p)
    :return:
    """

    ref_gene_info = {}
    for line in reference_info:
        col = line.split(',')
        vc, start, end, strand, chrom = col[2], col[3], col[4], col[5], col[6]
        ref_gene_info[vc] = {'start':start, 'end':end, 'strand':strand, 'chromosome':chrom}

    single_gene_sizes_list = []
    IGR_sizes= []
    for gene in genes_list:
        start = int(ref_gene_info[gene]['start'])
        end = int(ref_gene_info[gene]['end'])
        size = end - start
        if '_' not in gene:
            single_gene_sizes_list.append(size)
        if '_' in gene:
            IGR_sizes.append(size)

    reference_genome_length = int(reference_genome_length)
    #calculating all core genome stats
    total_core_genes_size = sum(single_gene_sizes_list)
    total_core_genome_percen = round((total_core_genes_size / reference_genome_length) * 100, 2)
    print('Total CG Sizes = ' + str(total_core_genes_size) + ' b.p (' + str(total_core_genome_percen) + '%)')

    ##calculating all IGR stats
    total_IGR_sizes = sum(IGR_sizes)
    total_IGR_percen = round((total_IGR_sizes / reference_genome_length) * 100, 2)
    print('Total IGR Sizes = ' + str(total_IGR_sizes) + ' b.p (' + str(total_IGR_percen) + '%)')

    ##calc the CG + IGR total
    total_IGR_and_CG = total_IGR_sizes + total_core_genes_size
    print('Total CG + IGR Size = ' + str(total_IGR_and_CG) + ' b.p')


    return total_core_genes_size, total_core_genome_percen, total_IGR_sizes, total_IGR_percen


calc_CG_bp_size(infile_list, reference_gff_info, reference_genome_length)