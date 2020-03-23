import argparse
import os
from time import sleep as sl
import pandas as pd
from random import random
pd.set_option('display.width', 400)
pd.set_option('display.max_columns', 15)
import random

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args

def missing_info_calc():

    #manually add to metadata file
    infile = open('/Users/liamcheneyy/Desktop/vcseventh_15/grapetree/non_grapetree_MGT9_allele_profiles.tsv','r').read().splitlines()
    results = {}
    for line in infile[1:]:
        no_allele_count = 0
        missing_info_count = 0
        col = line.split('\t')
        strain = col[0].split('.')[0]
        for cell in col[3:]:
            if '-' in cell:
                missing_info_count = missing_info_count + 1
            if cell == '0':
                no_allele_count = no_allele_count + 1
        results[strain] = [missing_info_count, no_allele_count]

    with open('/Users/liamcheneyy/Desktop/out.tsv','w') as out:
        out.write('Strain ID' + '\t' + 'Missing Info Allele' + '\t' + 'Failed Allele' + '\n')
        for key, value in results.items():
            out.write(key + '\t')
            for el in value:
                out.write(str(el) + '\t')
            out.write('\n')

def subset_df(mgt_want_level, df):

    # wont include strains with missing year, country, negative or missing alleles from MGT9 STDST or negative missing from seperate calculation
    # sub = df[(df['Year'].notnull()) & (df['Country'].notnull()) & (df['Negative or Missing Alleles'] == 0) & (~df['Country'].astype(str).str.contains('NaN') & df['MGT9 STDST'].astype(str).str.contains('\.0'))]

    # wont include strains with missing data or MGT9STDST
    sub = df[(df['Missing Info Allele'] <= 5) & (df['Failed Allele'] == 0)]

    sub = sub[["Strain ID", mgt_want_level,"Year", "Missing Info Allele", "Failed Allele"]]

    return sub

def filter_sts(ST_threshold, mgt_want_level, sub_df):

    #get frequency of each ST
    st_list = list(sub_df[mgt_want_level])
    st_freq_dict = {x:st_list.count(x) for x in st_list}

    # remove STs for small numbers of strains

    keep_dict = {}
    for key, value in st_freq_dict.items():
        if value >= ST_threshold:
            keep_dict[key] = value

    st_total_num = len(st_freq_dict.keys())
    chosen_st_num = len(keep_dict.keys())
    print("Total Number of STs: " + str(st_total_num))
    print("Number of Chosen STs: " + str(chosen_st_num))
    print("With a ST Threshold of " + str(ST_threshold) + ", number of missed STs is: " + str(st_total_num - chosen_st_num))
    print("Total Number of Selected Strains : " + str(sum(st_freq_dict.values())))
    return keep_dict

def seperate_sts(st_freq_dict,mgt_want_level,df):
    sep_st_dict = {}
    for item in st_freq_dict.keys():
        sub_df = df[df[mgt_want_level] == item]
        sub_df = sub_df.reset_index(drop=True)
        sep_st_dict[item] = sub_df

    return sep_st_dict

def organising_subdataframes(multiple_st_dict):

    for key, value in multiple_st_dict.items():
        value = value.sort_values('Year')
        value = value.reset_index(drop=True)
        multiple_st_dict[key] = value

    return multiple_st_dict

def selecting_strains(multiple_st_dict):

    keeping_list = []
    for key, value in multiple_st_dict.items():

        if value.shape[0] > 5:

            #select and append two strains with time
            sub_value = value[value['Year'].notnull()]

            if sub_value.shape[0] > 2:
                sub_value = sub_value.reset_index(drop=True)
                min_row_index = value['Year'].idxmin()
                max_row_index = value['Year'].idxmax()

                min_row = list(value.loc[min_row_index, :])
                max_row = list(value.loc[max_row_index, :])

                keeping_list.append(min_row)
                keeping_list.append(max_row)

            #select a random strain with no year metadata
            no_time_value = value[~value['Year'].notnull()]
            no_time_value = no_time_value.reset_index(drop=True)
            if no_time_value.shape[0] > 2:
                random_num = random.randint(0, no_time_value.shape[0] - 1)
                keep_row = list(no_time_value.loc[random_num, :])
                keeping_list.append(keep_row)

        #if ST smaller than 5 then randomly take sample
        else:
            random_num = random.randint(0, value.shape[0] - 1)
            keep_row = list(value.loc[random_num, :])
            keeping_list.append(keep_row)


    return keeping_list

    #organise using MGTST

def write_out(st_dict, out_directory, mgt_want_level, using_strains_df):
    if not os.path.exists(out_directory):
        os.mkdir(out_directory)

    for key, value in st_dict.items():
        with open(out_directory + '/' + mgt_want_level.split(' ')[0] + '_ST' + str(key) + '.txt','w') as out:
            for el in value:
                out.write(el + '.fasta' + '\n')

    using_strains_df.to_csv(out_directory + '/selected_strains.csv', sep=',', index=False)

def main():

    args = parseargs()
    df = pd.read_csv(df_in, sep='\t')

    ##get missing allele info
    # sub_df = missing_info_calc()

    #get wanted MGT leve, year and country data
    sub_df = subset_df(mgt_want_level, df)

    #selecting the STs
    keeping_st_dict = filter_sts(ST_threshold, mgt_want_level, sub_df)

    #create a dict with different subdataframes per MGT ST
    multiple_st_dict = seperate_sts(keeping_st_dict,mgt_want_level,sub_df)

    #sorting the subdataframes
    organised_multiple_st_dict = organising_subdataframes(multiple_st_dict)

    #selecting data
    keep_strains_list = selecting_strains(organised_multiple_st_dict)

    #write out strains lists to seperate files
    #write out sub_df to file not included
    # write_out(keep_strains_list, out_directory, mgt_want_level)

    for i in keep_strains_list:
        print(i[0])



if __name__ == '__main__':
    mgt_want_level  = "MGT6 ST"
    df_in = '/Users/liamcheneyy/Desktop/vcseventh_15/grapetree/all_metadata.txt'
    represent_df = '/Users/liamcheneyy/Desktop/representative/selected_strains.csv'
    out_directory = '/Users/liamcheneyy/Desktop/representative'
    ST_threshold = 2
    main()


###Idea to get strains from each of the major clades of the seventh pandemic tree.
##Used a desired MGT level.

#pre-calculation step.
#Get MGT9 allele allignment for all strains and calculate number of alleles with missing or negative info

##main calculation:
# choose STs assigned to > X strains.