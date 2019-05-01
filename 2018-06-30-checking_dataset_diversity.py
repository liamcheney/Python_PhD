import pandas as pd
from scipy import stats
import numpy as np
import random
import math
from time import clock
import openpyxl
import sys

df = pd.read_csv(sys.argv[1], low_memory=False)
outfile = sys.argv[2]
percentage = float(sys.argv[3])

df = pd.read_csv('/Users/liam/Desktop/96_gene_presence_absence_ref.csv', low_memory=False)
outfile = '/Users/liam/Desktop/output.csv'
percentage = 1

###isolate rows of DF which have genes found in 99% of strains
def df_isolate_random(df, percentage):
    df.insert(loc=0, column='Gene Count', value=df.notnull().sum(axis=1))
    strains_used = df.shape[1] - 1
    df = df.loc[df['Gene Count'] >= math.ceil(strains_used * percentage)]
    return df

def df_diversity_test(dataframe, out):
    results_df = pd.DataFrame()
    # for x in (i for j in (range(1,11,1), range(20,101,10)) for i in j):
    for x in range(1,11,1):
        result = []
        print(str(clock()) + '\t' + 'Running test at ' + str(x) + ' %.')
        for sample in range(1,5,1):
            strains_df = dataframe.iloc[:, 14:]
            percen = x / 100
            strains_df = strains_df.sample(frac = percen, axis=1)
            strains_df = df_isolate_random(strains_df, percentage)
            result.append(strains_df.shape[0])
        result.append(np.mean(result))
        result.append(stats.sem(result))
        results_df[str(percen) + '_' + str(strains_df.shape[1])] = result
    results_df.to_csv(out, index=False)
    return

df_diversity_test(df, outfile)


