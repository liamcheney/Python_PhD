from time import sleep as sl
import sys
import os
import math
import pandas as pd
import numpy as np
from time import clock
import glob
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import progressbar
import json
from io import StringIO
import csv

#reading in roary gene absence and presence
roary_path = "/Users/liamcheneyy/Desktop/roary/gene_presence_absence.csv"
gff_folder_in = "/Users/liamcheneyy/Desktop/roary/prokka_gffs/"
reference_accession = "GCA000006745"
info_dict_out = "/Users/liamcheneyy/Desktop/roary/info_dict.txt"
read_genome_info_dict = True
outfile_path = "/Users/liamcheneyy/Desktop/roary/"
isolate_percentage = 0.99

#dont need
write_out_roary_details = True
fix_roary_csv_val = False
write_out_fix_roary = False
read_in_fix_roary = False

##genome dictionary
def reading_or_creating_genomes_dict(genome_info_dict, info_dict_out):
    ###will read in previoulsy created dictionary or create genome info dict if needed
    if genome_info_dict == True:
        print(str(clock()) + '\t' + 'Reading in dictionary with all gff genome information.')
        with open(info_dict_out) as dict_file:
            info_dict = json.loads(dict_file.read())
            return info_dict
    else:
        print(str(clock()) + '\t' + 'Creating dictionary with all gff genome information.')
        sl(1)
        info_dict = creating_genome_info(gff_folder_in, reference_accession, info_dict_out)
        return info_dict
def creating_genome_info(gff_folder_in, reference_accesion, info_dict_out):
    number_of_genomes = len(list(glob.iglob(gff_folder_in + '/*.gff'))) + 1
    bar = progressbar.ProgressBar(maxval=number_of_genomes, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()
    info_dict = {}
    genome_count = 0
    for filename in glob.iglob(gff_folder_in + '/*.gff'):
        bar.update(genome_count)
        genome_count = genome_count + 1
        with open(filename, 'r') as file:
            strain_name = filename.split('/')[-1].rstrip('.gff')
            info_dict[strain_name] = {}
            for line in file:

                if ('CDS' in line and 'ID=' in line) or ('RNA' in line and 'ID=' in line):
                    col = line.split('\t')
                    beg_cord = col[3]
                    end_cord = col[4]
                    strand = col[6]
                    size = int(end_cord) - int(beg_cord)
                    gene_id = col[8].split(';')[0].split('_')[-1]
                    # if gene_id[-2] == '.':
                    #     gene_id = gene_id.split('.')[0]

                    for q in col[8].split(';'):
                        if "locus_tag" in q:
                            genbank = q.split('=')[-1].split('_')[-1]

                    if '_' not in col[0]:
                        contig = col[0]

                    elif '_' in col[0]:
                        if len(col[0].split('_')) > 2:
                            contig = col[0].split('_')[0] + '_' + col[0].split('_')[1]
                        else:
                            contig = col[0]

                    gene_id = gene_id.strip('=ID')

                    info_dict[strain_name][gene_id] = {'size': size, 'start': beg_cord, 'end': end_cord, 'contig': contig, 'locus_tag':genbank, 'strand':strand}
    bar.finish()

    #writing out the genome dictionary to file to save time
    print('\n')
    with open(info_dict_out, 'w') as dict_file:
        dict_file.write(json.dumps(info_dict))

    genome_info_dict =  True
    return info_dict
def gathering_core_gene_information(temp_file, info_dict):

    ###script will go over earch ortholog, first find average size of non paralogs, then find combined sizes of each paralog.
    print('\n')
    print(str(clock()) + '\t' + 'Gathering core gene information.')
    sl(1)
    percen_len = sum(1 for row in temp_file)
    bar = progressbar.ProgressBar(maxval=percen_len,
                                  widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()
    sl(1)
    group_count = 1

    for line in range(1, len(temp_file)):
    # for line in range(1, 6):

        ##creating progress var
        bar.update(group_count)
        group_count = group_count + 1

        # for each strain for the core gene
        for j in range(15, len(temp_file[line])):
            info = temp_file[line][j]
            # will leave empty strains
            if 'nan' in str(info):
                pass

            # handling cells which DONT have orthologs
            elif '\t' not in str(info):
                acc, contig, gene_id, beg_cord, end_cord, size, strand = handling_roary_single_annoations_strings(str(info), info_dict,reference_accession)

                temp_file[line][j] = (acc + '_' + str(contig) + '_' + str(gene_id) + '_' + str(beg_cord) + '_' + str(
                    end_cord) + '_' + str(size) + '(' + str(strand) + ')')

            # handling cells which have paralogs
            elif '\t' in str(info):
                frag_cell, frag_length_list = handling_roary_paralogs_annotations_strings(str(info), info_dict, reference_accession)
                temp_file[line][j] = frag_cell

    return temp_file

#fix roary naming
def fix_roary_csv(temp_file):

    # converting temp_file dataframe into list of lists for pythonic str functions etc
    temp_file = [temp_file.columns.values.tolist()] + temp_file.values.tolist()

    if fix_roary_csv_val:
        for line_num in range(1, len(temp_file), 1):
            for cell_num in range(15, len(temp_file[line_num]), 1):
                #fix .fasta in cells

                try:
                    if '.fasta' in temp_file[line_num][cell_num]:
                        temp_file[line_num][cell_num] = temp_file[line_num][cell_num].replace('.fasta','')

                except:
                    pass

                try:
                    #fix cds in cells
                    if 'cds-' in temp_file[line_num][cell_num]:
                        temp_file[line_num][cell_num] = temp_file[line_num][cell_num].replace('cds','')
                except:
                    pass

                try:
                # fix cds in cells
                    if '__' in temp_file[line_num][cell_num]:
                        temp_file[line_num][cell_num] = temp_file[line_num][cell_num].replace('__', '_')
                except:
                    pass

    if write_out_fix_roary:
        print("Writing out fixed Roary presence and absence.")
        with open(outfile_path + '/fixed_roary.csv', 'w') as out:
            for line in temp_file:
                for cell in line:
                    if 'nan' in str(cell):
                        out.write(',')
                    else:
                        out.write(str(cell) + ',')
                out.write('\n')

    return temp_file

#removing fragmented genomes with high paralogs
def isolate_ortho(roary_path, outfile_path, isolate_percentage):
    # function will isolate genes found in 99% of strains for each ortholog group
    print(str(clock()) + '\t' + 'Reading in ' + roary_path)
    df = pd.read_csv(roary_path, low_memory=False, index_col=False)
    print(str(clock()) + '\t' + 'Isolating orthologs with >99% frequency across all strains.')
    print('\n')

    strains_used = df.shape[1] - 13
    df.insert(loc = 6, column = 'Gene Count', value = df.iloc[:,13:].notnull().sum(axis=1))
    df = df.loc[df['Gene Count'] >= math.ceil(strains_used * isolate_percentage)]

    working_df = pd.DataFrame
    working_df = df.loc[df['Gene Count'] >= math.ceil(strains_used * isolate_percentage)]
    working_df.to_csv(outfile_path + '/filtered_gap.csv',sep=',')

    return working_df
def calculate_non_paralogous_core_genes(temp_file):
    # calculate the number of paralogous core genes for each genomes
    print(str(clock()) + '\t' + 'Calculating number of paralogous core genes per genome')
    genome_paralog_count_df = pd.DataFrame(columns=["Paralogous Core Gene Count"])
    for column in temp_file.columns.values[15:]:
        genome_paralog_count_df.loc[column] = temp_file[column].str.contains("\t").sum()

    #organising number of paralogs per genome dataframe in ascending order
    # genome_paralog_count_df = genome_paralog_count_df.set_index(list(genome_paralog_count_df)[0])
    genome_paralog_count_df = genome_paralog_count_df.sort_values(by="Paralogous Core Gene Count", ascending=True)

    percent_include = 100
    genomes_used_list = []
    non_para_core_genes_list = []
    paralogous_genomes_dict = {}

    #TODO
    for i in range(5,20,5):
        print(str(clock()) + '\t' + 'Calculating core genes after removing ' + str(i) + '% of highly paralogous genomes.')

        #with each iteration a list of genomes is created by slicing the dataframe using a range function.
        para_list_length = genome_paralog_count_df.shape[0]
        slice_length = int((percent_include / 100) * para_list_length)
        slice_df = genome_paralog_count_df[0:slice_length]
        sliced_genomes_list = slice_df.index.tolist()
        paralogous_genomes_dict[percent_include] = temp_file.columns.values.tolist()[0:15] + sliced_genomes_list

        #the sliced list is used to create a dataframe subset
        subset_roary_gap = temp_file[temp_file.columns.intersection(sliced_genomes_list)]

        # the number of non-paralogous core genes is recalculated based on the sliced list of genomes.
        total_non_para_core_gene_count = 0
        for index, row in subset_roary_gap.iterrows():
            para_core_genes_count = row.str.contains('\t').sum()
            non_para_core_genes_count = len(row) - para_core_genes_count
            if non_para_core_genes_count == len(row):
                total_non_para_core_gene_count = total_non_para_core_gene_count + 1

        genomes_used_list.append(percent_include)
        non_para_core_genes_list.append(total_non_para_core_gene_count)
        percent_include = 100 - i



    #TODO
    # ##creating a graph to show user the effects of removing genomes highest in paralogs on the non-paralogous core gene number
    # print('\n')
    # print(str(clock()) + '\t' + 'Presenting graph of non-paralogous core genes per genome.')
    # plt.close("all")
    # core_gene_data = {"Percentage of Used Genomes": genomes_used_list, "Number of Non-Paralogous Core Genes":non_para_core_genes_list}
    # graph_df = pd.DataFrame(data=core_gene_data)
    # x = graph_df[graph_df.columns[0]]
    # y = graph_df[graph_df.columns[1]]
    # plt.plot(x,y, '-o')
    # plt.xticks(np.arange(0,105,step=5))
    # plt.show()

    # #asking user input to choose the number of genomes to remain after removing genomes high in paralogs
    # plt.close("all")
    # print('\n')
    # # roary_to_keep_input = int(input("Enter percentage of strains to remain in analysis:"))
    # roary_to_keep_input = 100
    # temp_file = temp_file[temp_file.columns.intersection(paralogous_genomes_dict[roary_to_keep_input])]
    # with open('/Users/liamcheneyy/Desktop/pari/qc_paralogs_out.txt', 'w') as out:
    #     for i in paralogous_genomes_dict[roary_to_keep_input]:
    #         out.write(i + '\n')

    return temp_file
def merge(intervals):
   if not intervals:
       return []
   data = []
   for interval in intervals:
       data.append((interval[0], 0))
       data.append((interval[1], 1))
   data.sort()

   merged = []
   stack = [data[0]]
   for i in range(1, len(data)):
       d = data[i]
       if d[1] == 0:
           # this is a lower bound, push this onto the stack
           stack.append(d)
       elif d[1] == 1:
           if stack:
               start = stack.pop()
           if len(stack) == 0:
               # we have found our merged interval
               merged.append((start[0], d[0]))
   return merged

#handling Roary paralog and ortholog issues
def handling_roary_single_annoations_strings(j, info_dict, reference_accession):
    #will handle all the different roary annotations to extract gene_id, start, end and size
    if len(j.split('_')[0]) == 2 and j.split('_')[0].isalpha():
        acc = j.split('_')[0] + '_' + j.split('_')[1]
        gene_id = j.rstrip('"').split('_')[-1]

        contig = info_dict[acc][gene_id]['contig']
        beg_cord = info_dict[acc][gene_id]['start']
        end_cord = info_dict[acc][gene_id]['end']
        size = info_dict[acc][gene_id]['size']

        return acc, contig, gene_id, beg_cord, end_cord, size

    else:
        acc = j.split('_')[0]
        gene_id = j.rstrip('"').split('_')[-1]

        contig = info_dict[acc][gene_id]['contig']
        beg_cord = info_dict[acc][gene_id]['start']
        end_cord = info_dict[acc][gene_id]['end']
        size = info_dict[acc][gene_id]['size']
        strand = info_dict[acc][gene_id]['strand']

        return acc, contig, gene_id, beg_cord, end_cord, size, strand
def handling_roary_paralogs_annotations_strings(j, info_dict, reference_accession):
    frag_length_list = []
    frag_count = 0
    para_count = j.count('\t') + 1
    frag_list = []

    if len(j.split('_')[0]) == 2 and j.split('_')[0].isalpha():
        fragments = j.split('\t')
        for frag in range(para_count):
            acc = j.split('_')[0] + '_' + j.split('_')[1]
            gene_id = fragments[frag].rstrip('"').split('_')[-1]

            contig = info_dict[acc][gene_id]['contig']
            beg_cord = int(info_dict[acc][gene_id]['start'])
            end_cord = int(info_dict[acc][gene_id]['end'])
            size = info_dict[acc][gene_id]['size']
            strand = info_dict[acc][gene_id]['strand']


            frag_list.append(
                acc + '_' + contig + '_' + gene_id + '_' + str(beg_cord) + '_' + str(end_cord) + '_' + str(size) + '(' + str(strand) + ')')
            frag_length_list.append((beg_cord, end_cord))

    else:
        fragments = j.split('\t')
        for frag in range(para_count):
            acc = fragments[frag].split('_')[0]
            gene_id = fragments[frag].rstrip('"').split('_')[-1]
            # print(fragments[frag])

            contig = info_dict[acc][gene_id]['contig']
            beg_cord = int(info_dict[acc][gene_id]['start'])
            end_cord = int(info_dict[acc][gene_id]['end'])
            size = info_dict[acc][gene_id]['size']
            strand = info_dict[acc][gene_id]['strand']


            frag_list.append(
                acc + '_' + contig + '_' + gene_id + '_' + str(beg_cord) + '_' + str(end_cord) + '_' + str(size) + '(' + str(strand) + ')')
            frag_length_list.append((beg_cord, end_cord))

    frag_cell = '\t'.join(frag_list)
    return frag_cell, frag_length_list
def fasely_split_ortholgs_filter(reference_accession, keep_core_gene, temp_file, line, fail_reason):
    reference_core_gene = ""
    cell_in_size = 0
    true_para = 0
    false_para = 0
    for j in range(15, len(temp_file[line]), 1):
        j = str(temp_file[line][j])

        #find reference gene size if not fragmented
        if (reference_accession in j) and ('\t' not in j) and ('nan' not in j):
            reference_core_gene = j
            reference_size = int(j.split('_')[-1].split('(')[0])

            # find size of other cells in the same core group
            for l in temp_file[line][15:]:
                l = str(l)

                # if other cells are a single fragment
                if '\t' not in l and 'nan' not in l and l != '':
                    other_size = int(l.split('_')[-1].split('(')[0])

                    # compare size of reference cell against other cells
                    ref_lower_range = int(0.8 * reference_size)
                    ref_upper_range = int(1.2 * reference_size)
                    if ref_lower_range <= other_size <= ref_upper_range:
                        cell_in_size = cell_in_size + 1

                    # para_ref_lower_range = int(1.8 * reference_size)
                    # para_ref_upper_range = int(2.2 * reference_size)
                    # if para_ref_lower_range <= other_size <= para_ref_upper_range:
                    #     true_para = true_para + 1


                # if other cells are fragmented
                if '\t' in l and 'nan' not in l:
                    frag_sizes_list = []
                    fragments = l.split('\t')
                    for frag in fragments:
                        frag_size = int(frag.split('_')[-1].split('(')[0])
                        frag_sizes_list.append(frag_size)
                    total_frag_size = sum(frag_sizes_list)

                    # compare size of reference cell against other cells
                    ref_lower_range = int(0.8 * reference_size)
                    ref_upper_range = int(1.2 * reference_size)
                    if ref_lower_range <= total_frag_size <= ref_upper_range:
                        cell_in_size = cell_in_size + 1

                    # para_ref_lower_range = int(1.8 * reference_size)
                    # para_ref_upper_range = int(2.2 * reference_size)
                    # if para_ref_lower_range <= total_frag_size <= para_ref_upper_range:
                    #     true_para = true_para + 1

        elif reference_accession in j and '\t' in j and 'nan' not in j:
            fail_reason = "No intact reference gene."

        elif reference_accession in j and 'nan' in j:
            fail_reason = "No available reference gene."

    if cell_in_size <= (0.98 * len(temp_file[line][15:])):
        keep_core_gene = False
        fail_reason = "Failed, not >= 98% of strains with fragments within 20% of reference core gene."

    return keep_core_gene, reference_core_gene, fail_reason
def fasely_joined_orthologs_filter(reference_accession, split_core_gene, temp_file, line):
    reference_core_gene = ""
    for j in range(15, len(temp_file[line]), 1):
        j = str(temp_file[line][j])
        if reference_accession in j and '\t' in j and 'nan' not in j:
            reference_core_gene = j
            ##find combined length of fragmented cells
            frag_sizes_list = []
            fragments = j.split('\t')
            for frag in fragments:
                frag_size = int(frag.split('_')[-1].split('(')[0])
                frag_sizes_list.append(frag_size)
            total_frag_size = sum(frag_sizes_list)

            ##check if single gene is in line and if that gene is similar size to combined fragments
            single_gene_present = False
            for q in temp_file[line][15:]:
                q = str(q)
                if '\t' not in q and 'nan' not in q and q != '':
                    single_gene_size = int(q.split('_')[-1].split('(')[0])
                    lower_limit_ref = int(0.8 * total_frag_size)
                    upper_limit_ref = int(1.2 * total_frag_size)
                    if lower_limit_ref <= single_gene_size <= upper_limit_ref:
                        single_gene_present = True

            # if single gene present, then treat core gene as fasely joined paralogs
            if single_gene_present == True:
                fragments = j.split('\t')
                # comparing the sizes of fasely joined paralogs to sizes in other genome cells
                # taking reference fragment size
                both_frags_in_size = 0
                for frag in fragments:
                    frag_in_size = 0
                    frag_size = int(frag.split('_')[-1].split('(')[0])
                    # taking other cells fragment size
                    for other_cells in temp_file[line][15:]:
                        other_cells = str(other_cells)
                        other_cells_fragments = other_cells.split('\t')
                        for other_cells_frags in other_cells_fragments:
                            if 'nan' not in other_cells_frags and other_cells_frags != '':
                                other_cell_size = int(other_cells_frags.split('_')[-1].split('(')[0])

                                # creating upper and lower limits for variation in size between ref and other cell fragments
                                lower_limit_ref = int(0.8 * frag_size)
                                upper_limit_ref = int(1.2 * frag_size)

                                # comparing each reference fragment against the size of other cell fragments
                                if lower_limit_ref <= other_cell_size <= upper_limit_ref:
                                    frag_in_size = frag_in_size + 1

                        # checking if the fragment is found across >= 99% of the test genomes
                    if frag_in_size >= (0.99 * len(temp_file[line][15:])):
                        both_frags_in_size = both_frags_in_size + 1

                # check if should split core genes
                if both_frags_in_size == len(fragments):
                    split_core_gene = True
    return split_core_gene, reference_core_gene
def handling_roary_core_gene_ortholog_paralogs(temp_file, reference_accession):
    print('\n')
    print(str(clock()) + '\t' + 'Handling erroneous core genes.')

    core_gene_list = []
    excluded_list = []
    split_count = 0
    ##for each line of the Roary input
    print(str(clock()) + '\t' + 'Identifying fasely split orthologs and fasely joined paralogs.')
    for line in range(1, len(temp_file),1):
        keep_core_gene = True
        split_core_gene = False
        fail_reason = ""

        # print(temp_file[line])
        ##handle fasely split orthologs
        keep_core_gene, reference_core_gene_fso, fail_reason = fasely_split_ortholgs_filter(reference_accession, keep_core_gene, temp_file, line, fail_reason)

        ##handle fasely joined orthologs
        split_core_gene, reference_core_gene_fjp = fasely_joined_orthologs_filter(reference_accession, split_core_gene, temp_file, line)

        #add line to new_temp_file
        if keep_core_gene == True and split_core_gene == False:
            core_gene_list.append(reference_core_gene_fso)

        ##spliting fasely joined cells and appending
        if split_core_gene == True:
            split_count = split_count + 1

            single_gene = reference_core_gene_fjp.split('\t')
            for j in single_gene:
                core_gene_list.append(j)

        if keep_core_gene == False and split_core_gene == False:
            excluded_list.append(temp_file[line])


    return core_gene_list, excluded_list

#writing out core genes
def core_gene_out(core_gene_list, outfile_path, reference_accession, gff_folder_in, excluded_list):
    print('\n')
    print(str(clock()) + '\t' + 'Writing out list of core genes to: ' + outfile_path + '/core_genes.csv')
    with open(outfile_path + '/core_genes.csv','w') as outfile:
        outfile.write(','.join(["Accession", "CDS", "Locus Tag", "Start", "End", "Length"]))
        outfile.write('\n')
        for line in core_gene_list:
            col = line.split('_')
            acc = reference_accession
            cds = col[-4]
            start = col[-3]
            end = col[-2]
            length = col[-1]
            vc = converting_cds_to_vc(line, gff_folder_in, reference_accession)
            outfile.write(','.join([acc,cds,vc,start,end,length]))
            outfile.write('\n')

    with open(outfile_path + '/excluded_list.csv','w') as outfile:
        for line in excluded_list:
            for cell in line:
                outfile.write(str(cell) + ',')
            outfile.write('\n')
def converting_cds_to_vc(core_gene_in, gff_folder_in, reference_accession):
    vc = ""
    core_gene_cds = core_gene_in.split('_')[-4]
    for filename in glob.iglob(gff_folder_in + "/" + reference_accession + '*.gff'):
        with open(filename) as gff_in:
            for line in gff_in:
                if "CDS	" in line:
                    col = line.split("\t")
                    ref_cds = col[8].split(';')[0].split('_')[0][3:]

                    # find the line which has the core gene cds in the gff
                    if ref_cds == core_gene_cds:
                        parent = col[8].split(';')[1].split('=')[-1]

                        # find the parent line in gff which has VC
                        with open(filename) as gff_in:
                            for parent_line in gff_in:
                                if "gene	" in parent_line:
                                    col_1 = parent_line.split("\t")
                                    parent_gene = col_1[8].split(';')[0].split('_')[0][3:]

                                    # finding the parent gene in previous line
                                    if parent_gene == parent:
                                        vc = col_1[8].split(';')[-1].split('=')[-1].strip('\n')

    return vc

#master function
def master_handling_roary_paralog_problems(outfile_path, read_genome_info_dict, info_dict_out, reference_accession, roary_path, isolate_percentage, write_out_roary_details):

    # ##isolate core genes found in >=99% of the dataset genomes
    temp_file = isolate_ortho(roary_path, outfile_path, isolate_percentage)

    ##calculate the number of non paralogs core genes based on included genomes
    # temp_file = calculate_non_paralogous_core_genes(temp_file)

    ##fix roary naming
    temp_file = fix_roary_csv(temp_file)

    # ##will read in previoulsy created dictionary or create genome info dict if needed
    info_dict = reading_or_creating_genomes_dict(read_genome_info_dict, info_dict_out)

    ##gathering information core gene groups
    temp_file = gathering_core_gene_information(temp_file, info_dict)
    if write_out_roary_details:
        print()
        print("writing out roary with extra information.")
        with open(outfile_path + '/gap_with_details.csv','w') as out:
            for line in temp_file:
                for j in line:
                    out.write(str(j) + ',')
                out.write('\n')

    # ##handling Roary fasely split orthologs and fasely joined paralogs
    core_gene_list, excluded_list = handling_roary_core_gene_ortholog_paralogs(temp_file, reference_accession)

    ##will write out a list of core genes
    core_gene_out(core_gene_list, outfile_path, reference_accession, gff_folder_in, excluded_list)

master_handling_roary_paralog_problems(outfile_path, read_genome_info_dict, info_dict_out, reference_accession, roary_path, isolate_percentage, write_out_roary_details)

#TODO must change to do all analysis from 5,15 currently
# with open('/Users/liamcheneyy/Desktop/temp_file_in.csv','w') as outfile:
#     for line in temp_file:
#         for j in line:
#             outfile.write(str(j) + ',')
#         outfile.write('\n')