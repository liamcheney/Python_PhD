import pandas as pd
import math
from time import clock
from time import sleep

#structure of program
#1. read in lengths and organise based on paralog groups size (descending).
#2. create lists based containing increase sizs of strain names, begginning with strain with most paralog groups
#3. feed these lists into script to remove strains and find non-paralogous genes.
# overall will perform many analysis showing the chance in core genes when removing strains high in paralog groups.
print(str(clock()) + '\t' + 'creating lists of strains to remove in analysis.')
df_lengths = pd.read_csv('/Users/liam/Desktop/D3A_genome_length_with_paralog_groups.csv', header=0)
df_lengths.sort_values(by = 'Paralog Groups', ascending = False, inplace=True)

#creates list of 7th pandemic genomes
df_7th = df_lengths.loc[df_lengths['Pandemic'] == 'Seventh']
list_7th = df_7th['Accession'].tolist()

test_infile = open('/Users/liam/Desktop/seventh-strains.txt', 'r').read().splitlines()
test_list = [i for i in test_infile]
count = 0
for i in test_list:
    if i in list_7th:
        count += 1

#creates lists ranging the size of the all strains. smaller lists have strains with most paralog groups
accession =  df_lengths['Accession'].tolist()
cuttoff_analysis = {}
for x in range(0, len(accession), int(0.05*len(accession))):
    cuttoff_analysis[x] = accession[:x]

## calculate the average number of paralogs per strains for each analysis
para_sum_list = []
used_list = []
length_and_para_dict = dict(zip(df_lengths['Accession'], df_lengths['Paralog Groups']))
for acc_list in cuttoff_analysis.values():
    reduced_list = sorted(set(used_list)^set(acc_list))
    para_sum = 0
    for j in reduced_list:
        para_sum = para_sum + length_and_para_dict[j]
        used_list.append(j)
    para_sum_list.append(para_sum)

# provide input strains to be removed (remove their columns) and script will remove and output.
print(str(clock()) + '\t' + 'loading gene absence and presence file and removing strains with most paralogs.')
df = pd.read_csv('/Users/liam/Desktop/gene_presence_absence.csv', index_col=0, low_memory=False)
graph_out = open('/Users/liam/Desktop/D2C_graph.csv', 'w')
graph_out.write('Core Gene Amount' + ',' + 'Strains Used in Analysis' + ',' + 'Percentage of Seventh Strains' + '\n')

para_list_count = 0
for j in cuttoff_analysis.values():
    para_list_count = para_list_count + 1
    count_7th = 0
    j = list(set(j))
    print(str(clock()) + '\t' + 'begginning analysis removing ' + str(len(j)) + ' strains')
    sub_df = df.drop(j, axis=1)

    #a collection of pandas fucntions. 1st will create a column counting how many genes a strain has. 2nd will keep all lines if a gene is present in 99 or more strains.
    print(str(clock()) + '\t' + 'isolating orthologs with >99% frequency across all strains.')
    strains_used = sub_df.shape[1] - 15
    sub_df.insert(loc = 6, column = 'Gene Count', value = sub_df.iloc[:,14:].notnull().sum(axis=1))
    kept_lines = sub_df.loc[sub_df['Gene Count'] >= math.ceil(strains_used * 0.99)].to_csv().splitlines()

    # calculates the amount of remaining strains in the seventh pandemic
    for n in j:
        if n in list_7th:
            count_7th = count_7th + 1
    count_7th = int((count_7th / len(j)) * 100)

    ##will select only lines which have no paralogs and write them to an outfile
    outfile = open('/Users/liam/Desktop/D2C_missing_' + str(len(j)) + '_strains_no_paralogs_presence_absence.csv', 'w')
    print(str(clock()) + '\t' + 'writing out orthologs with no paralogs present')
    CG_count = 0
    for i in kept_lines:
        gene_count = 0
        paralog_count = 0
        group_name = i.split(',')[0]
        strains_cols = i.split(',')[14:]
        for k in strains_cols:
            if k.count('\t') == 0:
                gene_count = gene_count + 1
            if k.count('\t') > 0:
                paralog_count = paralog_count + 1

        if paralog_count == 0:
            # outfile.write(i + '\n')
            CG_count = CG_count + 1

    graph_out.write(str(CG_count) + ',' + str(len(strains_cols)) + ',' + str(para_sum_list[para_list_count]) + '\n')
    # outfile.close()
    print('\n')

# sub process for removing strains from a gap which are paralogous, however will not remove the paralog lines afterwards.
#used mainly for looking at annotation issues causing high paralog numbers

import pandas
from time import clock
import math

df = pandas.read_csv('/Users/liam/Desktop/D3A_i92_GAP.csv')
infile = open('/Users/liam/Desktop/strains.txt','r').read().splitlines()

drop_list = [x for x in infile]

print(str(clock()) + '\t' + 'begginning analysis removing ' + str(len(infile)) + ' strains')
sub_df = df.drop(drop_list, axis=1)

#a collection of pandas fucntions. 1st will create a column counting how many genes a strain has. 2nd will keep all lines if a gene is present in 99 or more strains.
print(str(clock()) + '\t' + 'isolating orthologs with >99% frequency across all strains.')
strains_used = sub_df.shape[1] - 14
sub_df.insert(loc = 6, column = 'Gene Count', value = sub_df.iloc[:,14:].notnull().sum(axis=1))
kept_lines = sub_df.loc[sub_df['Gene Count'] >= math.ceil(strains_used * 0.99)].to_csv('/Users/liam/Desktop/D3A_i92_missing175strains_withparalogs.csv')
