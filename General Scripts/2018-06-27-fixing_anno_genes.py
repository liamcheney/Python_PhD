import glob
import pandas as pd
from time import sleep as sl

infile = open('/Users/liam/Desktop/99percen_96i_para_ana_gap.csv', 'r').read().splitlines()
outfile = open('/Users/liam/Desktop/99percen_96i_para_gaps_ana_gap.csv', 'w')

outfile.write('1bp' + ',' + '2-5bp' + ',' + '6-10bp' + ',' + '11-20bp' + ',' + '21-30bp' + ',' + '31-40bp' + ',' + '41-50bp' + ',' + '51-100bp' + ',' + '>101bp' + ',' + 'Cell Count' + ',' + 'Paralog Cell Count' + ',' + 'More Than Two Paralogs' + ',' + 'Overlapping' + '\n')
for j in infile[1:]:
    firstbp = 0
    secondbp = 0
    thirdbp = 0
    fourthbp = 0
    fifthbp = 0
    sixbp = 0
    sevenbp = 0
    eightbp = 0
    ninebp = 0
    para_cell_count = 0
    cell_count = 0
    more_than_two_para = 0
    overlap = 0

    strains = j.split(',')[14:]
    for k in strains:
        cell_count = cell_count + 1
        if '\t' in k:
            para_cell_count = para_cell_count + 1
            cord_list = []
            para_count = k.count('\t') + 1
            fragment = k.split('\t')
            for frag in range(para_count):
                acc = fragment[frag].split('_')[0]
                gene_id = fragment[frag].split('_')[3]
                node = fragment[frag].split('_')[1] + '_' +fragment[frag].split('_')[2]
                start = fragment[frag].split('_')[-3]
                end = fragment[frag].split('_')[-2]
                cord_list.append(start)
                cord_list.append(end)

            if len(cord_list) == 4:
                calc = int(cord_list[2]) - int(cord_list[1])
                if calc < 0:
                    overlap = overlap + 1

                if calc == 1:
                    firstbp = firstbp + 1
                if calc >= 2  and calc <= 5:
                    secondbp = secondbp + 1
                if calc > 6 and calc <= 10:
                    thirdbp = thirdbp + 1
                if calc > 11 and calc <= 20:
                    fourthbp = fourthbp + 1
                if calc > 21 and calc <= 30:
                    fifthbp = fifthbp + 1
                if calc > 31 and calc <= 40:
                    sixbp = sixbp + 1
                if calc > 41 and calc <= 50:
                    sevenbp = sevenbp + 1
                if calc > 51 and calc <= 100:
                    eightbp = eightbp + 1
                if calc > 101:
                    ninebp = ninebp + 1

            if len(cord_list) > 4:
                more_than_two_para  = more_than_two_para + 1



    outfile.write(str(firstbp) + ',' + str(secondbp) + ',' + str(thirdbp) + ',' + str(fourthbp) + ',' + str(fifthbp) + ',' + str(sixbp) + ',' + str(sevenbp) + ',' + str(eightbp) + ',' + str(ninebp) + ',' + str(cell_count) + ',' + str(para_cell_count) + ',' + str(more_than_two_para) + ',' + str(overlap) + '\n')