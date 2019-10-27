import argparse
from time import sleep as sl
import pandas as pd

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args

def choose_level_least_STs(keep_dict):

    save_list = []
    for key, value in keep_dict.items():
        if len(value) != 0:
            save_list.append([key, len(value)])

    save_list.sort(key=lambda x: x[1])
    chosen_level = save_list[0][0]
    keep_dict = keep_dict[save_list[0][0]]

    return keep_dict, chosen_level

def remove_non_specific_st(st_save_dict, min_strains_per_st, percentage_of_contaminating_strains):

    # remove any STs which are not specific
    keep_dict = {}
    for mgt_level in st_save_dict:
        keep_dict[mgt_level] = {}

        for ST in st_save_dict[mgt_level]:

            total_false_true_calls = st_save_dict[mgt_level][ST]['true'] + st_save_dict[mgt_level][ST]['false']
            contamiating_strains = int(percentage_of_contaminating_strains / 100 * total_false_true_calls)

            # if an ST covers both true and false
            if st_save_dict[mgt_level][ST]['true'] >= 1 and st_save_dict[mgt_level][ST]['false'] >= 1:

                # if ST shares a small amount of true and false strains, then remain, not removed.
                if (st_save_dict[mgt_level][ST]['true'] <= contamiating_strains or st_save_dict[mgt_level][ST]['false'] <= contamiating_strains):
                    keep_dict[mgt_level][ST] = st_save_dict[mgt_level][ST]

                    # if ST has true and false within accepted issues, filter out smaller STs
                    if st_save_dict[mgt_level][ST]['true'] >= min_strains_per_st or st_save_dict[mgt_level][ST][
                        'false'] >= min_strains_per_st:
                        keep_dict[mgt_level][ST] = st_save_dict[mgt_level][ST]

            # if ST has no overlap of TRUE and FALSE, filter out smaller STs
            elif st_save_dict[mgt_level][ST]['true'] >= min_strains_per_st or st_save_dict[mgt_level][ST]['false'] >= min_strains_per_st:
                keep_dict[mgt_level][ST] = st_save_dict[mgt_level][ST]

    return keep_dict

def filtering_level(st_save_dict, strains_number, min_strains_per_st, percentage_of_contaminating_strains):

    # remove any STs which are not specific
    keep_dict = remove_non_specific_st(st_save_dict, min_strains_per_st, percentage_of_contaminating_strains)

    # choose the level with least STs
    keep_dict, chosen_level = choose_level_least_STs(keep_dict)

    # choose the STs describing user defined % the majority of clade
    final_st_dict = {}
    true_strain_count = 0
    true_st_list = []
    false_strain_count = 0
    false_st_list = []

    for key,value in keep_dict.items():

        if value['true'] > value['false']:

            true = value['true']
            true_st_list.append(key)
            true_strain_count = true_strain_count + true

        elif value['false'] > value['true']:

            false = value['false']
            false_st_list.append(key)
            false_strain_count = false_strain_count + false

    true_percentage_st_reps = round((true_strain_count / strains_number * 100),2)
    false_percentage_st_reps = round((false_strain_count / strains_number * 100),2)

    print(chosen_level, 'True', true_st_list, true_strain_count, true_percentage_st_reps)
    print(chosen_level, 'False', false_st_list, false_strain_count, false_percentage_st_reps)

    return final_st_dict

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
    start_col = 191
    end_col = 192
    min_strains_per_st = 10
    percentage_of_contaminating_strains = 5

    #read in metadata
    df = pd.read_csv('/Users/liamcheneyy/Desktop/vcseventh_21/grapetree/seventh/metadata.csv', index_col='ID', low_memory=False)
    strains_number = df.shape[0]
    #columns of attributes to take
    want_attributes = list(df.columns.values[start_col:end_col])

    #will return MGT level and best ST to describe group
    save_dict = {}
    for element in want_attributes:
        print(element)

        # calculate the STs for each attribute
        st_save_dict = sts_per_attributes(element, df)

        # select the best level to define attribute
        filtered_st_n_level = filtering_level(st_save_dict, strains_number, min_strains_per_st, percentage_of_contaminating_strains)

if __name__ == '__main__':
    main()



#Logic
#to find the STs which best describe attributes

#challenges
#take the wanted gene,
#for each level calculate the

