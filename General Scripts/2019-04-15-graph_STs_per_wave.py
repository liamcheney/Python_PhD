from time import sleep as sl
import pandas as pd

infile_path = "/Users/liamcheneyy/Desktop/select_Metadata_grapetree.txt"
list_of_levels = ['MGT2 ST', 'MGT3 ST', 'MGT4 ST']
out_path = "/Users/liamcheneyy/Desktop/select_wave_STs.tsv"

def read_in(infile_path):

    df = pd.read_csv(infile_path, sep='\t')

    wave1_df = df[df['Wave'] == 1.0]
    wave2_df = df[df['Wave'] == 2.0]
    wave3_df = df[df['Wave'] == 3.0]

    return df, wave1_df, wave2_df, wave3_df

def make_st_waves(wave1_df, wave2_df, wave3_df):
    wave1_freq_dict = {}

    for col_num in range(0, len(list_of_levels), 1):

        col_head = "MGT" + str(col_num + 2) + " ST" ###made since there is no MGT1
        wave1_freq_dict[col_head] = list(wave1_df[list_of_levels[col_num]].value_counts())

    wave2_freq_dict = {}

    for col_num in range(0, len(list_of_levels), 1):
        col_head = "MGT" + str(col_num + 2) + " ST"  ###made since there is no MGT1
        wave2_freq_dict[col_head] = list(wave2_df[list_of_levels[col_num]].value_counts())

    wave3_freq_dict = {}

    for col_num in range(0, len(list_of_levels), 1):
        col_head = "MGT" + str(col_num + 2) + " ST"  ###made since there is no MGT1
        wave3_freq_dict[col_head] = list(wave3_df[list_of_levels[col_num]].value_counts())

    return wave1_freq_dict, wave2_freq_dict, wave3_freq_dict


def main(infile_path, list_of_levels, out_path):

    #read in file and create a dict with all info
    df, wave1_df, wave2_df, wave3_df = read_in(infile_path)

    #frequencies for each dict
    wave1_freq_dict, wave2_freq_dict, wave3_freq_dict = make_st_waves(wave1_df, wave2_df, wave3_df)

    with open(out_path,'w') as out:
        out.write('Wave 1' + '\n')
        for key in wave1_freq_dict:
            for el in wave1_freq_dict[key]:
                out.write(str(el) + '\t')
            out.write('\n')

        out.write('Wave 2' + '\n')
        for key in wave2_freq_dict:
            for el1 in wave2_freq_dict[key]:
                out.write(str(el1) + '\t')
            out.write('\n')

        out.write('Wave 3' + '\n')
        for key in wave3_freq_dict:
            for el2 in wave3_freq_dict[key]:
                out.write(str(el2) + '\t')
            out.write('\n')


main(infile_path, list_of_levels, out_path)