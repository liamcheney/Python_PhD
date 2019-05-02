from time import sleep as sl
import glob
import pandas as pd
import sys
import numpy as np
import json
import time
#
# #creates dictionary with all gene information for each strain
# info_dict = {}
# for filename in glob.iglob('/Users/liam/Desktop/2018-04-02-Gffs_D2C/*.gff'):
#     with open(filename, 'r') as file:
#         strain_name = filename.split('/')[-1].rstrip('.gff')
#         info_dict[strain_name] = {}
#         print(strain_name)
#         if 'GCA_000006745' in filename:
#             for line in file:
#                 if 'CDS' in line and 'ID=' in line:
#                     col = line.split('\t')
#                     contig = col[0]
#                     gene_id = col[8].split(';')[0].split('_')[0][3:]
#
#                     beg_cord = col[3]
#                     end_cord = col[4]
#                     size = int(end_cord) - int(beg_cord)
#
#                     info_dict[strain_name][gene_id] = {'size': size, 'start': beg_cord, 'end': end_cord, 'contig': contig}
#         elif 'GC' in filename:
#             for line in file:
#                 if 'CDS' in line and 'ID=' in line:
#                     col = line.split('\t')
#                     contig = col[0]
#                     gene_id = col[8].split(';')[0][-5:]
#
#                     beg_cord = col[3]
#                     end_cord = col[4]
#                     size = int(end_cord) - int(beg_cord)
#
#                     info_dict[strain_name][gene_id] = {'size': size, 'start': beg_cord, 'end': end_cord, 'contig': contig}
#
#         elif 'ER' in filename:
#             for line in file:
#                 if 'CDS' in line and 'ID=' in line:
#                     col = line.split('\t')
#                     contig = col[0].split('_')[1] + '_' + col[0].split('_')[2]
#                     gene_id = col[8].split(';')[0][-5:]
#
#                     beg_cord = col[3]
#                     end_cord = col[4]
#                     size = int(end_cord) - int(beg_cord)
#
#                     info_dict[strain_name][gene_id] = {'size':size, 'start':beg_cord, 'end':end_cord, 'contig':contig}
#
# with open('/Users/liam/Desktop/D2C_geneid_dict_with_cord.txt', 'w') as dict_file:
#     dict_file.write(json.dumps(info_dict))

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

info_dict = {}
with open('/Users/liam/Desktop/D2C_geneid_dict_with_cord.txt') as dict_file:
    info_dict = json.loads(dict_file.read())

# script will go over earch ortholog, first find average size of non paralogs, then find combined sizes of each paralog.
infile = open('/Users/liam/Desktop/D2C_i96_missing345strains_withparalogs.csv', 'r').read().splitlines()
outfile = open('/Users/liam/Desktop/D2C_i96_99percen_para_ana.csv','w')

outfile.write(infile[0] + ',' + ' Single Gene Mean' + ',' + 'Single Gene Average' + ',' + 'Small False Paralogs' + ',' + 'False Paralogs' + ',' + 'Indeterminable Paralogs' + ',' + 'True Paralogs' + ',' + 'Large True Paralogs' + ',' + 'Single Genes' + ',' + 'Single Gene Small' + '\n')

for i in infile[1:]:
    frag_count_list = []
    single_length_count_list = []

    true_para = 0
    false_para = 0
    large_true_para = 0
    inde_para = 0
    small_false_para = 0
    single_gene_length = 0

    col = i.split(',')
    strains = col[15:]
    for q in col[0:15]:
        outfile.write(q + ',')
    for j in strains:
        #will write out empty strains
        if 'ERR' not in j and 'GC' not in j and 'cds' not in j:
            outfile.write(',')

        ##as input files name format differs, have to run similar code three times to get information from both ERR and GCA groups.
        elif '\t' not in j and j != '':
            if 'cds' in j:
                col = j.split('_')
                acc = 'GCA_000006745'
                gene_id = col[0].lstrip('"')

                beg_cord = info_dict[acc][gene_id]['start']
                end_cord = info_dict[acc][gene_id]['end']
                size = info_dict[acc][gene_id]['size']
                single_length_count_list.append(size)
                outfile.write(gene_id + '_' + str(beg_cord) + '_' + str(end_cord) + '_' + str(size) + ',')

            if 'ERR' in j:
                acc = j.split('_')[0].lstrip('"').rstrip('"')
                gene_id = j.rstrip('"')[-8:-3]

                contig = info_dict[acc][gene_id]['contig']
                beg_cord = info_dict[acc][gene_id]['start']
                end_cord = info_dict[acc][gene_id]['end']
                size = info_dict[acc][gene_id]['size']
                single_length_count_list.append(size)
                outfile.write(acc + '_' + contig + '_' + gene_id + '_' + beg_cord + '_' + end_cord + '_' + str(size) + ',')

            if 'GC' in j:
                acc = j.split('_')[0:2]
                if '.' in acc[1]:
                    acc = 'GCA_' + acc[1][:-2] # need to convert all to GCA for my files
                else:
                    acc = 'GCA_' + acc[1]
                gene_id = j.rstrip('"')[-8:-3]

                contig = info_dict[acc][gene_id]['contig']
                beg_cord = info_dict[acc][gene_id]['start']
                end_cord = info_dict[acc][gene_id]['end']
                size = info_dict[acc][gene_id]['size']

                single_length_count_list.append(size)
                outfile.write(acc + '_' + contig + '_' + gene_id + '_' + beg_cord + '_' + end_cord + '_' + str(size) + ',')

        elif '\t' in j and j != '':
            frag_length_list = []
            frag_count = 0
            para_count = j.count('\t') + 1
            frag_list = []

            if 'cds' in j:
                fragments = j.split('\t')
                for frag in range(para_count):
                    col = fragments[frag].split('_')
                    acc = 'GCA_000006745'
                    gene_id = col[0].lstrip('"')
                    beg_cord = int(info_dict[acc][gene_id]['start'])
                    end_cord = int(info_dict[acc][gene_id]['end'])
                    size = info_dict[acc][gene_id]['size']

                    frag_list.append(gene_id + '_' + str(beg_cord) + '_' + str(end_cord) + '_' + str(size))
                    frag_length_list.append((beg_cord, end_cord))

            if 'ERR' in j:
                fragments = j.split('\t')
                for frag in range(para_count):
                    acc = fragments[frag].split('_')[0].lstrip('"')
                    gene_id = fragments[frag].split('_')[-1][0:-3].rstrip('(')

                    contig = info_dict[acc][gene_id]['contig']
                    beg_cord = int(info_dict[acc][gene_id]['start'])
                    end_cord = int(info_dict[acc][gene_id]['end'])
                    size = info_dict[acc][gene_id]['size']

                    frag_list.append(acc + '_' + contig + '_' + gene_id + '_' + str(beg_cord) + '_' + str(end_cord) + '_' + str(size))
                    frag_length_list.append((beg_cord,end_cord))

            if 'GC' in j:
                fragments = j.split('\t')
                for frag in range(para_count):
                    acc = j.split('_')[0:2]
                    if '.' in acc[1]:
                        acc = 'GCA_' + acc[1][:-2] # need to convert all to GCA for my files
                    else:
                        acc = 'GCA_' + acc[1]

                    gene_id = fragments[frag].split('_')[-1][0:-3].rstrip('(')
                    contig = info_dict[acc][gene_id]['contig']
                    beg_cord = int(info_dict[acc][gene_id]['start'])
                    end_cord = int(info_dict[acc][gene_id]['end'])
                    size = info_dict[acc][gene_id]['size']

                    frag_list.append(acc + '_' + contig + '_' + gene_id + '_' + str(beg_cord) + '_' + str(end_cord) + '_' + str(size))
                    frag_length_list.append((beg_cord, end_cord))

            frag_cell = '\t'.join(frag_list)
            outfile.write(frag_cell + ',')

            if len(frag_length_list) >= 2:
                cords_list = merge(frag_length_list)
                total_size = sum([(size[1] - size[0]) for size in cords_list])
                frag_count_list.append(total_size)

    ###for calculating the types of paralogs inside an ortholog group
    mean = np.mean(single_length_count_list)
    median = np.median(single_length_count_list)
    for length in frag_count_list:
        if length <= (0.8 * mean):
            small_false_para = small_false_para + 1
        if (0.8 * mean) <= length <= (1.2 * mean):
            false_para = false_para + 1
        if (1.2 * mean) <= length <= (1.8 * mean):
            inde_para = inde_para + 1
        if (1.8 * mean) <= length <= (2.2 * mean):
            true_para = true_para + 1
        if length >= (2.2 * mean):
            large_true_para = large_true_para + 1

    for lengths in single_length_count_list:
        if lengths <= (0.8 * mean):
            single_gene_length = single_gene_length + 1

    outfile.write(str(mean) + ',' + str(median) + ',' + str(small_false_para) + ',' + str(false_para) + ',' + str(inde_para) + ',' + str(true_para) + ',' + str(large_true_para) + ',' + str(len(single_length_count_list)) + ',' + str(single_gene_length) + '\n')
outfile.close()