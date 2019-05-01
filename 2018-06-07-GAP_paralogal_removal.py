import pandas as pd
import math
from time import clock
from time import sleep
import glob

graph_list = []
for filename in glob.iglob('/Users/liam/Desktop/IGR_gap/*'):
    print(str(clock()) + '\t' + 'Reading in ' + filename)
    df = pd.read_csv(filename, low_memory=False)
    blast_percentage = filename.split('/')[-1][0:2]
    print(str(clock()) + '\t' + 'Isolating orthologs with >99% frequency across all strains.')
    strains_used = df.shape[1] - 15
    df.insert(loc = 6, column = 'Gene Count', value = df.iloc[:,14:].notnull().sum(axis=1))
    df= df.loc[df['Gene Count'] >= math.ceil(strains_used * 0.99)]
    kept_lines = df.loc[df['Gene Count'] >= math.ceil(strains_used * 0.99)].to_csv().splitlines()

    # outfile = open('/Users/liamcheneyy/Desktop/IGR_no_paralogs/' + blast_percentage + '_no_paralogs_presence_absence.csv', 'w')

    print(str(clock()) + '\t' + 'Writing out orthologs with no paralogs present')
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

    graph_list.append(str(CG_count) + ',' +  str(len(strains_cols)) + ',' + blast_percentage + '\n')
    # outfile.close()
    print('\n')

graph_out = open('/Users/liam/Desktop/IGR_gap/D1A_noparalogs_graph.csv', 'w')
graph_out.write('Core Gene Amount' + ',' + 'Strains Used in Analysis' + '\n')

for x in graph_list:
    graph_out.write(x)
