import argparse
from time import sleep as sl
import pandas as pd

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args

def remove_non_specific_st(st_save_dict, min_strains_per_st, percentage_of_contaminating_strains):

    # remove any STs which are not specific
    keep_dict = {}
    for mgt_level in st_save_dict:
        keep_dict[mgt_level] = {}

        for ST in st_save_dict[mgt_level]:

            total_false_true_calls = st_save_dict[mgt_level][ST]['true'] + st_save_dict[mgt_level][ST]['false']
            contaminating_strains = int(percentage_of_contaminating_strains / 100 * total_false_true_calls)

            #if an ST has more false than true cant be used
            if st_save_dict[mgt_level][ST]['false'] > st_save_dict[mgt_level][ST]['true']:
                pass

            # if an ST covers both true and false
            elif st_save_dict[mgt_level][ST]['true'] >= 1 and st_save_dict[mgt_level][ST]['false'] >= 1:

                # if ST shares a small amount of false strains, then remain, not removed.
                # print(mgt_level, str(ST), str(total_false_true_calls), contaminating_strains)
                if st_save_dict[mgt_level][ST]['false'] <= contaminating_strains:
                    keep_dict[mgt_level][ST] = st_save_dict[mgt_level][ST]

                    # if ST has true and false within accepted issues, filter out smaller STs
                    if st_save_dict[mgt_level][ST]['true'] >= min_strains_per_st:
                        keep_dict[mgt_level][ST] = st_save_dict[mgt_level][ST]

            # if ST has no overlap of TRUE and FALSE, filter out smaller STs
            elif st_save_dict[mgt_level][ST]['true'] >= min_strains_per_st:
                keep_dict[mgt_level][ST] = st_save_dict[mgt_level][ST]

    return keep_dict

def find_best_level_and_ST(element, strains_with_attribute_list, st_save_dict):

    save_dict = {}

    wanted_number = len(strains_with_attribute_list)

    #find number of true per level
    save_list = []
    for key, value in st_save_dict.items():
        print(key, value)
        true_count = 0
        st_num_count = 0
        for i in value:
            true_count = true_count + st_save_dict[key][i]['true']
            st_num_count = st_num_count + 1
        save_list.append([key, true_count, st_num_count])

    #organise list by MOST number of trues and then LEAST ST number
    save_list.sort(key= lambda x: (x[1], -x[2]), reverse=True)
    want_level = save_list[0][0]
    print(save_list)

    percentage_got = round(save_list[0][1] / wanted_number * 100,1)

    save_dict['level'] = want_level
    save_dict['percentage described'] = percentage_got
    save_dict['total number'] = wanted_number
    save_dict['STs'] = list(st_save_dict[want_level].keys())

    print("Chosen level :", want_level)
    print("Percentage Described :", str(percentage_got))


    return save_dict

def sts_per_attributes(element, df, min_strains_per_st):

    #create sub dataframe of MGT ST levels and attribute
    sub_df_cols = [element] + [x for x in df if 'MGT' in x and 'CC' not in x and 'MGT9' not in x]
    sub_df = df[sub_df_cols]
    mgt_col_list = [x for x in sub_df if "MGT2" in x or "MGT3" in x or "MGT4" in x]

    st_save_dict = {}
    for item in mgt_col_list:
        # print("Calculating STs per scheme for : ", item)
        want_levels = [sub_df_cols[0], item]
        sub_level_df = sub_df[want_levels]
        st_freq_dict = sub_level_df[item].value_counts().to_dict()

        st_save_dict[item] = {}
        for key, value in st_freq_dict.items():
            true_sub = sub_level_df[(sub_level_df[element] == True) & (sub_level_df[item] == key)]
            true_count = true_sub.shape[0]
            true_strains_list = list(true_sub.index)

            false_sub = sub_level_df[(sub_level_df[element] == False) & (sub_level_df[item] == key)]
            false_count = false_sub.shape[0]
            false_strains_list = list(false_sub.index)

            if true_count > 0:
                st_save_dict[item][key] = {'true':true_count, 'false':false_count, 'true_strains_list':true_strains_list, 'false_strains_list':false_strains_list}

    return st_save_dict

def strains_with_attribute(df, element):

    sub_df = df[df[element] == True]
    sub_df = list(sub_df[element].index)

    return sub_df

def main():
    args = parseargs()

    #variables
    start_col = 32
    end_col = 33
    min_strains_per_st = 10
    percentage_of_contaminating_strains = 10

    #read in metadata
    df = pd.read_csv('/Users/liamcheneyy/Desktop/vcseventh_22/grapetree/seventh/MGT_isolate_data.txt', index_col='ID', low_memory=False, delimiter='\t')

    #columns of attributes to take
    want_attributes = list(df.columns.values[start_col:end_col])

    #will return MGT level and best ST to describe group

    for element in want_attributes:
        print("Element :", element)

        #calculate number of strains with atrribute
        strains_with_attribute_list = strains_with_attribute(df, element)

        #calculate the STs for each attribute
        st_save_dict = sts_per_attributes(element, df, min_strains_per_st)

        #remove non specific STs
        st_save_dict = remove_non_specific_st(st_save_dict, min_strains_per_st, percentage_of_contaminating_strains)

        #select the best level to define attribute
        find_best_level_and_ST_dict = find_best_level_and_ST(element, strains_with_attribute_list, st_save_dict)

if __name__ == '__main__':
    main()



#Logic
#to find the STs which best describe attributes

#1. get a list of all strains TRUE for gene
#2. then make dict, for each MGT level, get the frequency of STs, the for each ST count the number of TRUE and FALSE for gene
#3. quality filter: remove small STs, remove ST which has too many false and true strains, need just TRUE within threshold
#4. organise the dict, and count
#
#
#

#challenges
#DONE: get list of strains with attribute
#DONE: got pandas calc for true and false each ST
#DONE: filter out non-specific STs

#compare the ST for each level against my list, take the ST which describes the most: problem: what if multiple STs?

#just add total number of trues for each level,