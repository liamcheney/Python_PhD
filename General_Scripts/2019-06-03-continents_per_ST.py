import pandas as pd
from time import sleep as sl
import collections
import csv

def read_in(infile_path):

    df = pd.read_csv(infile_path, sep='\t')
    return df

def ST_in_continent(df, STs_series, mgt_want_lvl):
    result_dict = {}

    #go over each ST for wanted MGT lvl
    for row, label in STs_series.iteritems():
        result_dict["ST" + str(row)] = {}

        #save the contients for each ST into a dict
        sub_df = df[df[mgt_want_lvl] == row]
        temp_list = sub_df['Country'].tolist()

        #dict to count
        counts = collections.Counter(temp_list)
        counts.most_common()
        counts = dict(counts)
        result_dict["ST" + str(row)] = counts

    return result_dict

def write_out(result_dict, outfile_path):


    with open(outfile_path,'w') as out:
        for key in result_dict:
            print(key)
            out.write(key + '\n')
            for key, value in result_dict[key].items():
                out.write(str(key) + ',' + str(value) + ',')
            out.write('\n')



def main(infile_path,outfile_path):

    mgt_want_lvl = 'MGT7 CC mergemin'

    #read in df
    df = read_in(infile_path)

    #STs for wanted mgt lvl
    STs_series = df[mgt_want_lvl].value_counts()

    #find how many continents each ST is found in
    result_dict = ST_in_continent(df, STs_series, mgt_want_lvl)

    #write out
    write_out(result_dict, outfile_path)

infile_path = "/Users/liamcheneyy/Desktop/metadata.tsv"
outfile_path = "/Users/liamcheneyy/Desktop/MGT7_CCMERGIN_continents_per_ST.csv"

main(infile_path, outfile_path)
