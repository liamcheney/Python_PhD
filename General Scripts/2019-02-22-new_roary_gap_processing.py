from time import sleep as sl
import sys
import math
import pandas as pd
from time import clock
import glob
import matplotlib.pyplot as plt

infile = sys.argv[1]

#function will isolate genes found in 99% of strains for each ortholog group
def isolate_ortho(gap_in):

    for filename in glob.iglob(gap_in):
        print(str(clock()) + '\t' + 'Reading in ' + filename)

        df = pd.read_csv(filename, low_memory=False)
        print(str(clock()) + '\t' + 'Isolating orthologs with >99% frequency across all strains.')

        strains_used = df.shape[1] - 15
        df.insert(loc = 6, column = 'Gene Count', value = df.iloc[:,14:].notnull().sum(axis=1))
        df= df.loc[df['Gene Count'] >= math.ceil(strains_used * 0.99)]
        working_df = df.loc[df['Gene Count'] >= math.ceil(strains_used * 0.99)]

        #calculate the number of paralogous core genes for each genomes
        print(str(clock()) + '\t' + 'Calculating number of paralogous core genes per genome')
        genome_paralog_count_df = pd.DataFrame(columns=["Paralogous Core Gene Count"])
        for column in working_df.columns.values[15:]:
            genome_paralog_count_df.loc[column] = working_df[column].str.contains("\t").sum()

        plt.close("all")
        genome_paralog_count_df.plot.hist(bins=25)
        string_name = int(input("cut off: "))
        plt.xlabel("Paralogous Core Gene Count")
        plt.legend().remove()
        plt.show()
        return working_df

isolate_ortho(infile)
