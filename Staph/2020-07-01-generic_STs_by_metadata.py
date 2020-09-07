#yeah ST above a certain minimum size that is >80% found in one country

import argparse
from time import sleep as sl
import multiprocessing as mp
import pandas as pd

def df_setup(infile_path,meta, assoc_meta):

    #read in metadata
    df = pd.read_csv(infile_path, index_col='Strain', low_memory=False, sep='\t')

    #select MGT columns and metadata
    col_list = []
    for column in df:
        if ('MGT' in column and 'CC' not in column):
            col_list.append(column)

    sub_df = df[col_list].astype(int)

    # df[list(df.columns)[0:8]] = df[list(df.columns)[0:8]].astype(int)
    if assoc_meta != '':
        want_col = [meta] + assoc_meta
    else:
        want_col = meta

    for el in want_col:
        sub_df[el] = df[el]
    print(sub_df.shape)
    return sub_df

def get_meta_variations(df):
    save_list = []
    meta_variants = df[meta].unique().tolist()
    for x in meta_variants:
        if not isinstance(x, float):
            save_list.append(x)
    return save_list

def relate_ST_to_meta(df, meta,assoc_meta, ST, level, ST_dict, max_contamination, min_st_size, done_list, save_dict):
    ST_size = ST_dict[ST]

    if ST_size > min_st_size:
        #create a sub dir for ST, find the meta dist for that ST on base leve
        sub_df = df[df[level] == ST]
        ST_size = sub_df.shape[0]
        meta_freq = sub_df[meta].value_counts().to_dict()
        largest_meta = list(meta_freq.keys())[0]
        largest_meta_val = meta_freq[largest_meta]
        min_meta_percen = 100 - max_contamination

        #if ST not contaminated
        meta_covered = round((largest_meta_val / ST_size * 100),2)
        if  meta_covered > min_meta_percen:
            MGT_ST = level + '_' + str(ST)
            if MGT_ST not in done_list:
                print(f'{level} ST{ST} {largest_meta_val} {meta_covered} {largest_meta}', sep='\t')

                want_col = [meta] + assoc_meta
                result_df = sub_df[want_col]
                result_df = result_df[result_df[meta] == largest_meta]
                save_dict[MGT_ST] = result_df
                done_list.append(MGT_ST)

        #need to break down ST using higher level
        else:
            next_level = 'MGT' + str(int(level[-1]) + 1)
            next_ST_dict = df[next_level].value_counts().to_dict()
            for next_ST in next_ST_dict:
                relate_ST_to_meta(df, meta,assoc_meta, next_ST, next_level, next_ST_dict, max_contamination, min_st_size, done_list, save_dict)

def iterator(df, meta, assoc_meta, max_contamination, min_st_size):
    # gets the lowest level, tries that, then further divides up
    print(f'MGT_level ST Num_Covered Percen_Covered Major_Country', sep='\t')
    seed_level = df.columns.tolist()[0]
    seed_ST_dict = df[seed_level].value_counts().to_dict()
    done_list = []
    save_dict = {}
    for seed_ST in seed_ST_dict:
        relate_ST_to_meta(df, meta,assoc_meta, seed_ST, seed_level, seed_ST_dict, max_contamination, min_st_size, done_list, save_dict)

    return save_dict

def out_ST_data(ST_df_dict, meta, assoc_meta):
    # print(f'ID MGTST Level ST Country Year', sep='\t')
    saves_list = []
    for key, value in ST_df_dict.items():
        for index,row in value.iterrows():
            year = list(row.values)[-1]
            country= list(row.values)[-2]
            ST = list(row.values)[0]
            level = key.split('_')[0]
            print(index, key, level, ST, country, year, sep='\t')


if __name__ == '__main__':
    #variables
    meta = 'Phage Genotype'
    # assoc_meta = ['']
    assoc_meta = ''
    infile_path = '/Users/liamcheneyy/Desktop/newgenotyping/genotyping_w_meta.txt'
    outfile_path = '/Users/liamcheneyy/Desktop/genotypes_by_ST.txt'
    min_st_size = 5
    max_contamination = 30
    get_ST_out = False

    ##set up the data frame
    df = df_setup(infile_path, meta, assoc_meta)
    df[meta] = df[meta].str.replace(' ','_')

    #related STs with Metadata
    ST_df_dict = iterator(df, meta, assoc_meta, max_contamination, min_st_size)

    # output isolates from STs
    if get_ST_out:
        out_ST_data(ST_df_dict, meta, assoc_meta)
