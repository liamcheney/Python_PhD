from time import sleep as sl
import glob
import pandas as pd

cds_input = open('/Users/liam/Desktop/D1A_cds_list','r').readlines()

for filename in glob.iglob('/Users/liam/Desktop/D1A_all_gab/40*csv'):
    fname = "blast precentage " + filename.replace('/Users/liam/Desktop/D1A_all_gab/','')[0:2]
    df = pd.read_csv(filename, index_col=0)
    ref_col = df['GCF_000006745']
    # for line in df.iterrows():
        

