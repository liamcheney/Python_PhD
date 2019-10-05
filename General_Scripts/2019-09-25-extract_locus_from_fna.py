from time import sleep as sl
from Bio import SeqIO
import sys

input_list = list(open(
    '/Users/liamcheneyy/Desktop/roary/core_genes_list.txt',
    'r').read().splitlines())
reference_in = open(
    '/Users/liamcheneyy/Desktop/roary/prokka_gffs/GCA000006745.gff',
    'r').read().splitlines()
reference_fasta = SeqIO.parse('/Users/liamcheneyy/Desktop/ref/GCA_000006745.fna',
                              'fasta')
output_dir = '/Users/liamcheneyy/Desktop/'


def reference_information_dict_creator(reference_in, table_or_gff):
    info_dict = {}

    if table_or_gff == "table":
        for i in reference_in:
            col = i.split('\t')
            if 'gene' in col[0]:
                locus_tag = col[16]
                start = col[7]
                end = col[8]
                chr = col[5]
                strand = col[9]
                length = int(end) - int(start)
                info_dict[locus_tag] = {'start': start, 'end': end, 'strand': strand, 'chromosome': chr, 'length': length}
    elif table_or_gff == "gff":
        for line in reference_in:
            if 'ID=' in line:
                col = line.split('\t')
                start = col[3]
                end = col[4]
                chr = col[0]
                strand = col[6]
                length = int(end) - int(start)
                locus_tag = int(col[8].split(';')[0].split('_')[-1])
                info_dict[locus_tag] = {'start': start, 'end': end, 'strand': strand, 'chromosome': chr, 'length': length}

    else:
        print('Error: Input annotation requires: "gff" or "table"')
        sys.exit()


    print(info_dict.keys())
    sl(1)
    return info_dict


def reading_input(input_list):
    save_list = []
    for line in input_list:
        if '_' in line:
            col = line.split('_')
            first = col[0]
            second = col[1]
            save_list.append([first, second])
        else:
            save_list.append([line])
    return save_list


def input_positions(positions_dict, input_save):
    save_dict = {}
    for element in input_save:
        if len(element) == 2:
            beggining_igr = element[0]
            end_igr = element[1]
            igr_start = int(positions_dict[beggining_igr]['end'])
            igr_end = int(positions_dict[end_igr]['start'])
            igr_name = beggining_igr + '_' + end_igr
            igr_length = igr_end - igr_start
            save_dict[igr_name] = {'start': igr_start, 'end': igr_end, 'length': igr_length}
        if len(element) == 1:
            locus = element[0]
            start = int(positions_dict[locus]['start'])
            end = int(positions_dict[locus]['end'])
            length = end - start
            save_dict[locus] = {'start': start, 'end': end, 'length': length}
    return save_dict


def extracting_sequences(input_positions_dict, reference_fasta):
    save_dict = {}
    for record in reference_fasta:
        chromosome = str(record.seq)
        for key, value in input_positions_dict.items():
            locus = key
            start = int(value['start']) - 1
            end = int(value['end']) - 1
            target_sequence = chromosome[start:end]
            save_dict[key] = target_sequence
    return save_dict


def write_out(input_positions_dict, sequences_dict, output_dir):
    print()
    # writing out fasta
    with open(output_dir + '/locus_regions.fasta', 'w') as out:
        for key, value in sequences_dict.items():
            out.write('>' + key + ':1' + '\n')
            out.write(value + '\n')
    # writing out details
    with open(output_dir + '/locus_positions.csv', 'w') as out:
        out.write('Locus_tag' + ',' + 'Start' + ',' + 'End' + ',' + 'Length' + '\n')
        for key, value in input_positions_dict.items():
            start = value['start']
            end = value['end']
            length = value['length']
            out.write(key + ',' + str(start) + ',' + str(end) + ',' + str(length) + '\n')


def main(reference_in, input_list, reference_fasta, output_dir):
    table_or_gff = "gff" #"table"

    # reading in gff information from FEATURE_TABLE
    positions_dict = reference_information_dict_creator(reference_in, table_or_gff)

    # read in input_list
    input_save = reading_input(input_list)
    # create positions from input dict
    input_positions_dict = input_positions(positions_dict, input_save)
    # create dictionary of input sequences
    sequences_dict = extracting_sequences(input_positions_dict, reference_fasta)
    # writing out
    write_out(input_positions_dict, sequences_dict, output_dir)


if __name__ == '__main__':
    main(reference_in, input_list, reference_fasta, output_dir)