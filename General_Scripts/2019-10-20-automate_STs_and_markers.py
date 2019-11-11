import argparse
from time import sleep as sl
import pandas as pd

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args

def sts_per_attributes(element, sub_df, min_strains_per_st, percen_contam_strains, final_save):
    save_list = []
    strains_list = []
    strains_described = 0
    total_true_strains = sub_df[sub_df[element] == True].shape[0]
    cols_list = sub_df.columns.values[0:-1]
    for column in cols_list:

        st_freq = sub_df[column].value_counts()
        st_freq_dict = st_freq[st_freq >= min_strains_per_st].to_dict()

        #calculate the False and True strains per ST of level
        for ST, COUNT in st_freq_dict.items():
            true_sub = sub_df[(sub_df[element] == True) & (sub_df[column] == ST)]
            true_count = true_sub.shape[0]

            false_sub = sub_df[(sub_df[element] == False) & (sub_df[column] == ST)]
            false_count = false_sub.shape[0]

            #if ST IS NOT contaminated then add to dict
            strains_with_st = COUNT
            if false_count <= (percen_contam_strains / 100 * strains_with_st):

                save_list.append(column + ' ST' + str(ST))
                strains_described = strains_described + true_count
                strains_list.append(list(true_sub.index))

                #remove any non comtaminated STs from sub_df
                sub_df = sub_df[sub_df[column] != ST]

            #if ST IS contaminated, then pass
            if false_count >= (percen_contam_strains / 100 * strains_with_st):
                pass
                # print(column, ST, COUNT, true_count, false_count)

    strains_list_final = [i for x in strains_list for i in x]
    percen_described = round((strains_described / total_true_strains * 100),2)
    num_of_st = len(save_list)
    final_save[element] = [percen_described] + [num_of_st] + save_list

    return final_save,strains_list_final


def main():
    args = parseargs()

    #variables
    start_col = 36
    end_col = 37
    min_strains_per_st = 10
    percen_contam_strains = 10

    #read in metadata
    df = pd.read_csv('/Users/liamcheneyy/Desktop/seventh_MGT_isolate_data.txt', index_col='ID', low_memory=False, delimiter='\t')

    #columns of attributes to take
    want_attributes = list(df.columns.values[start_col:end_col])

    ##will return MGT level and best ST to describe group
    final_save = {}
    strains_list = []

    ##calculting non-overlapping STs
    print("Element, Percentage, MGT STs".format('\t'))
    for element in want_attributes:

        ##saving dict
        final_save[element] = {}

        ##create sub_df of element and STs
        col_list = list(df.columns)[0:8] + [element]
        sub_df = df[col_list]

        #TODO edit for both T and F

        #calculate the STs for each attribute
        st_save_dict,strains_list = sts_per_attributes(element, sub_df, min_strains_per_st, percen_contam_strains, final_save)

    # for key, value in final_save.items():
    #     print(key, *value, sep='\t')

        ##creating graphing data

        #select strains from the previous ST
        #select strains with year metadata
        return_df = df[(df.index.isin(strains_list)) & (df['Year'].notnull())]
        # return_df = return_df[list(df.columns)[0:8] + ['Year']]

        #get the times for each ST
        mgt_sts_and_years = {}
        for key, value in st_save_dict.items():
            mgt_sts_and_years[key] = {}
            ST_list = value[2:]
            for el in ST_list:
                level = el.split(' ')[0]
                ST = el.split(' ')[-1].strip('ST')

                years_df = return_df[return_df[level] == float(ST)]
                years_df = years_df[[level] + ['Year']]
                years_freq = years_df['Year'].value_counts().to_dict()

                mgt_sts_and_years[key][level + ' ST' + str(ST)] = years_freq

        #graph out changes
        for key, value in mgt_sts_and_years.items():
            print(key,value)

        # print(return_df['Year'].value_counts())
        # print(sorted(list(set(return_df['Year']))))
if __name__ == '__main__':
    main()



#Logic
#to find the STs which best describe attributes

#1. get a list of all strains TRUE and FALSE for gene
#2. then make dict, for each MGT level, get the frequency of STs, the for each ST count the number of TRUE and FALSE for gene
#3. quality filter: remove small STs, remove ST which has too many false and true strains, need just TRUE within threshold
#4. organise the dict, and count
#
#
#


#old algorythm to find a single best level for an attribute
#didnt work since never clean definitions.


#def remove_non_specific_st(st_save_dict, min_strains_per_st, percentage_of_contaminating_strains):
    #
    # # remove any STs which are not specific
    # keep_dict = {}
    # for mgt_level in st_save_dict:
    #     keep_dict[mgt_level] = {}
    #
    #     for ST in st_save_dict[mgt_level]:
    #
    #         total_false_true_calls = st_save_dict[mgt_level][ST]['true'] + st_save_dict[mgt_level][ST]['false']
    #         contaminating_strains = int(percentage_of_contaminating_strains / 100 * total_false_true_calls)
    #
    #         #if an ST has more false than true cant be used
    #         if st_save_dict[mgt_level][ST]['false'] > st_save_dict[mgt_level][ST]['true']:
    #             pass
    #
    #         # if an ST covers both true and false
    #         elif st_save_dict[mgt_level][ST]['true'] >= 1 and st_save_dict[mgt_level][ST]['false'] >= 1:
    #
    #             # if ST shares a small amount of false strains, then remain, not removed.
    #             # print(mgt_level, str(ST), str(total_false_true_calls), contaminating_strains)
    #             if st_save_dict[mgt_level][ST]['false'] <= contaminating_strains:
    #                 keep_dict[mgt_level][ST] = st_save_dict[mgt_level][ST]
    #
    #                 # if ST has true and false within accepted issues, filter out smaller STs
    #                 if st_save_dict[mgt_level][ST]['true'] >= min_strains_per_st:
    #                     keep_dict[mgt_level][ST] = st_save_dict[mgt_level][ST]
    #
    #         # if ST has no overlap of TRUE and FALSE, filter out smaller STs
    #         elif st_save_dict[mgt_level][ST]['true'] >= min_strains_per_st:
    #             keep_dict[mgt_level][ST] = st_save_dict[mgt_level][ST]
    #
    # return keep_dict

#def find_best_level_and_ST(element, strains_with_attribute_list, st_save_dict):

    # save_dict = {}
    #
    # wanted_number = len(strains_with_attribute_list)
    #
    # #find number of true per level
    # save_list = []
    # for key, value in st_save_dict.items():
    #     print(key, value)
    #     true_count = 0
    #     st_num_count = 0
    #     for i in value:
    #         true_count = true_count + st_save_dict[key][i]['true']
    #         st_num_count = st_num_count + 1
    #     save_list.append([key, true_count, st_num_count])
    #
    # #organise list by MOST number of trues and then LEAST ST number
    # save_list.sort(key= lambda x: (x[1], -x[2]), reverse=True)
    # want_level = save_list[0][0]
    # print(save_list)
    #
    # percentage_got = round(save_list[0][1] / wanted_number * 100,1)
    #
    # save_dict['level'] = want_level
    # save_dict['percentage described'] = percentage_got
    # save_dict['total number'] = wanted_number
    # save_dict['STs'] = list(st_save_dict[want_level].keys())
    #
    # print("Chosen level :", want_level)
    # print("Percentage Described :", str(percentage_got))
    #
    #
    # return save_dict