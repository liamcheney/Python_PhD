import pandas as pd
import math
from time import clock
from time import sleep
import glob
from functions import isolate_ortho_cd



#function will isolate genes found in 99% of strains for each ortholog group
def isolate_ortho_cd(gap_in, gap_out, graph_out, percentage=99):
    graph = open(graph_out, 'a')
    graph.write('Blast Percentage' + ',' + 'Core Gene Number' + '\n')

    for filename in glob.iglob(gap_in):
        print(str(clock()) + '\t' + 'Reading in ' + filename)

        df = pd.read_csv(filename, low_memory=False)
        print(str(clock()) + '\t' + 'Isolating orthologs with >99% frequency across all strains.')

        strains_used = df.shape[1] - 15
        df.insert(loc = 6, column = 'Gene Count', value = df.iloc[:,14:].notnull().sum(axis=1))
        df= df.loc[df['Gene Count'] >= math.ceil(strains_used * percentage)]
        kept_lines = df.loc[df['Gene Count'] >= math.ceil(strains_used * percentage)].to_csv(gap_out, index=False)
        graph.write(str(percentage) + ',' + str(kept_lines.shape[0]) + '\n')


