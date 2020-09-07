import argparse
from time import sleep as sl
from pandas import ExcelWriter
import pandas as pd
import matplotlib
def remove_small_sts(result_dict, min_st_size, percen_out):

    # convert to list
    for k, v in result_dict.items():
        save_list = []
        collapse_num = 0
        y = list(v)
        total_st = sum(y)
        for x in y:
            if x >= min_st_size:
                if percen_out:
                    q = round((x / total_st * 100),2)
                    save_list.append(q)
                else:
                    save_list.append(x)
            else:
                collapse_num = collapse_num + x

        if percen_out:
            collapse_num = round((collapse_num / total_st * 100), 2)
            save_list.append(collapse_num)
        else:
            save_list.append(collapse_num)

        result_dict[k] = save_list

    return result_dict
def convert_to_percen(result_dict):

    save_dict = {}
    for k, v in result_dict.items():
        total_ST = v.shape[0]
        x = v.apply(lambda x: ((x / total_ST) * 100))
        x = x.round(2)
        save_dict[k] = x

    return save_dict
def get_sts_per_scheme(df):

    result_dict = {}
    for col in df:
        value_count = df[col].value_counts()
        result_dict[col] = value_count
        print(col, len(list(set(value_count))))

    return result_dict
def read_in(input):

    df = pd.read_csv(input, sep='\t')
    keep_col_list = [x for x in df.columns.values if 'MGT' in x and 'DST' not in x and 'CC' not in x]
    df = df[df.columns.intersection(keep_col_list)]
    # df = df[(~df['MGT2'].str.contains('None')) & (~df['MGT3'].str.contains('None')) & (~df['MGT4'].str.contains('None')) & (~df['MGT5'].str.contains('None')) & (~df['MGT6'].str.contains('None')) & (~df['MGT7'].str.contains('None')) & (~df['MGT8'].str.contains('None'))]
    return df

def main():
    #I/O
    input = '/Users/liamcheneyy/Desktop/MGT_isolate_data.txt'
    output = '/Users/liamcheneyy/Desktop/STs_per_scheme.txt'

    #convert to percentages
    percen_out = False

    #min_st_size
    min_st_size = 5

    #read in
    df = read_in(input)

    #get STs per scheme
    result_dict = get_sts_per_scheme(df)

    #collapse smaller STs
    result_dict = remove_small_sts(result_dict, min_st_size, percen_out)
    for k,v in result_dict.items():
        print(k, *v, sep='\t')
        # print(sum(v))

if __name__ == '__main__':
    main()
