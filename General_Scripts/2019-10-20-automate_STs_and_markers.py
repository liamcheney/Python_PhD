import argparse
from time import sleep as sl
import pandas as pd

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args

def filtering_level(st_save_dict, strains_number, describe_clade_percen):

    # remove any STs which are not specific
    keep_dict = {}
    for mgt_level in st_save_dict:
        keep_dict[mgt_level] = {}

        for ST in st_save_dict[mgt_level]:
            if st_save_dict[mgt_level][ST]['true'] > 0 and st_save_dict[mgt_level][ST]['false'] > 0:
                pass
            else:
                keep_dict[mgt_level][ST] = st_save_dict[mgt_level][ST]

    # choose the level with least STs
    save_list = []
    for key, value in keep_dict.items():
        save_list.append([key, len(value)])

    save_list.sort(key=lambda x: x[1])
    keep_dict = keep_dict[save_list[0][0]]

    # choose the STs describing 95% the majority of clade
    final_st_dict = {}
    want_strain_number = int(describe_clade_percen * strains_number / 100)
    strain_count = 0
    st_list = []
    for key,value in keep_dict.items():
        if strain_count < want_strain_number:
            true = value['true']
            st_list.append(key)
            strain_count = strain_count + true

    percentage_st_reps = round((strain_count / strains_number * 100),2)
    print(st_list, strain_count, percentage_st_reps)

def sts_per_attributes(element, df):

    print("Calculating STs per scheme")
    #create sub dataframe of MGT ST levels and attribute
    sub_df_cols = [element] + [x for x in df if 'MGT' in x and 'CC' not in x and 'MGT9' not in x]
    sub_df = df[sub_df_cols]
    mgt_col_list = [x for x in sub_df if "MGT2" in x or "MGT3" in x]

    st_save_dict = {}
    for item in mgt_col_list:
        want_levels = [sub_df_cols[0], item]
        sub_level_df = sub_df[want_levels]
        st_freq_dict = sub_level_df[item].value_counts().to_dict()

        st_save_dict[item] = {}
        for key, value in st_freq_dict.items():
            true_count = sub_level_df[(sub_level_df[element] ==  True) & (sub_level_df[item] == key)].shape[0]
            false_count = sub_level_df[(sub_level_df[element] ==  False) & (sub_level_df[item] == key)].shape[0]

            st_save_dict[item][key] = {'true':true_count, 'false':false_count}

    return st_save_dict

def main():
    args = parseargs()

    #variables
    start_col = 28
    end_col = 29
    include_CCs = False
    describe_clade_percen = 95

    #read in metadata
    df = pd.read_csv('/Users/liamcheneyy/Desktop/metadata.csv', index_col='ID', low_memory=False)
    strains_number = df.shape[0]
    #columns of attributes to take
    want_attributes = list(df.columns.values[start_col:end_col])

    #will return MGT level and best ST to describe group
    save_dict = {}
    for element in want_attributes:

        # calculate the STs for each attribute
        st_save_dict = sts_per_attributes(element, df)

        # select the best level to define attribute
        filtered_st_n_level = filtering_level(st_save_dict, strains_number, describe_clade_percen)

if __name__ == '__main__':
    main()

#Logic
#to find the STs which best describe attributes

#challenges
#2. the script must choose the best MGT level
    #find the level which defines the least amount of STs, with the largest number.
#3. script must choose an ST unique to a certain attribute




