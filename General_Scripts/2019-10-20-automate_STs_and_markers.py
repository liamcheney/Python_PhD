import argparse
from time import sleep as sl
import pandas as pd

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args
def df_setup(infile_path, start_col,end_col):

    #read in metadata
    df = pd.read_csv(infile_path, index_col='ID', low_memory=False, delimiter='\t')
    df.replace("None.None", "0", inplace=True)

    df[list(df.columns)[0:8]] = df[list(df.columns)[0:8]].astype(int)

    #columns of attributes to take
    want_attributes = list(df.columns.values[start_col:end_col])

    return df, want_attributes
def sts_per_attributes(element, sub_df, min_strains_per_st, percen_contam_strains,precen_inconsistent_strains,relate_results_to, df):
    mgt_st_list = []
    strains_list = []
    number_TRUE = 0
    total_strains = sub_df[sub_df[element] == True]
    total_strains_num = total_strains.shape[0]
    cols_list = sub_df.columns.values[0:-1]
    st_freq_dict = {}

    for column in cols_list:

        sub_st_freq = sub_df[column].value_counts()
        sub_st_freq = sub_st_freq[sub_st_freq >= min_strains_per_st].to_dict()

        #calculate the False and True strains per ST of level
        for ST, COUNT in sub_st_freq.items():
            sub = sub_df[(sub_df[element] == True) & (sub_df[column] == ST)]
            True_count = sub.shape[0]

            false_sub = sub_df[(sub_df[element] == False) & (sub_df[column] == ST)]
            false_count = false_sub.shape[0]

            #if ST IS NOT contaminated then add to dict
            strains_with_st = COUNT

            #each ST reset hierarch flag
            if false_count <= (percen_contam_strains / 100 * strains_with_st):

                #check the ST has not hierach inconc
                hierach_check = check_hierarch_incon(column, ST, precen_inconsistent_strains,sub)
                if hierach_check == "Pass":

                    mgt_st_list.append(column + ' ST' + str(ST))
                    number_TRUE = number_TRUE + True_count
                    strains_list.append(list(sub.index))

                    st_freq_dict[column + ' ST' + str(ST)] = {'st_count':True_count, 'genome_list':list(sub.index)}

                    #remove any non comtaminated STs from sub_df
                    sub_df = sub_df[sub_df[column] != ST]

            #if ST IS contaminated, then pass
            if false_count >= (percen_contam_strains / 100 * strains_with_st):
                pass

    percen_with_described = round((number_TRUE / total_strains_num * 100),2)
    num_of_st = len(mgt_st_list)
    all_strains = df.shape[0]
    percen_all_strains = round((number_TRUE / all_strains * 100),2)
    non_described_strains = total_strains_num - number_TRUE
    final_save = {'st_num':num_of_st, 'all_strains_in_dataset':all_strains, 'percen_all_strains':percen_all_strains, 'number_TRUE': number_TRUE, 'percen_with_described':percen_with_described, 'non_described_strains':non_described_strains, 'mgt_sts':st_freq_dict}

    strains_list = [i for x in strains_list for i in x]

    return_df = df[df.index.isin(strains_list)]
    not_st_strains_list = list(total_strains[~total_strains.index.isin(strains_list)].index)
    not_st_strains_df = df[df.index.isin(not_st_strains_list)]

    return final_save, return_df, not_st_strains_df
def check_hierarch_incon(column, ST, precen_inconsistent_strains, sub):

    #for a ST, will get the STs from previous level that describe ST.
    #if many previous level STs describe current ST then is hierachael incocnsistency

    if column == "MGT2":
        return "Pass"
    else:

        previous_level = "MGT" + str(int(column.replace("MGT",""))-1)
        prev_level_sts = sub[previous_level]
        strains_num = sub.shape[0]
        allowed_precen_inconsistent_strains = round(( precen_inconsistent_strains / 100 ) * strains_num, 0)
        prev_level_sts_freqs = prev_level_sts.value_counts().to_dict()

        #count the number of inconsistencies
        inconsis_count = 0
        for i in list(prev_level_sts_freqs.values())[1:]:
            inconsis_count = inconsis_count + i

        #if inconsistencies less than allowed number, then return
        if inconsis_count <= allowed_precen_inconsistent_strains:
            return "Pass"
def st_freq_and_waves(return_df, st_save_dict, element, not_st_strains_df, min_st_for_figure):
    # relate STs freq to waves
    wave_list = [1, 2, 3]
    waves_st_dict = {}
    waves_st_failed = {}

    #breaking down STs per Wave
    for el in wave_list:
        waves_st_dict[el] = {}
        for MGTST in st_save_dict['mgt_sts']:
            level = MGTST.split(' ')[0]
            ST = int(MGTST.split(' ')[1].strip('ST'))
            sub_df = return_df[(return_df['Wave'] == el) & (return_df[level] == ST)]

            if sub_df.shape[0] != 0:
                waves_st_dict[el][MGTST] = sub_df.shape[0]

    total_classified = 0
    total_unclassified = 0
    st_num = 0
    print(element, st_save_dict)
    #printing out STs per wave
    for k, v in waves_st_dict.items():
        for i, x in v.items():
            if x >= min_st_for_figure:
                ST = i.split(' ')[-1]
                gene = element.split('_')[0]
                wave = "Wave " + str(k)
                freq = x
                level = i.split(' ')[0]
                total_classified = total_classified + freq
                st_num = st_num + 1
                print(gene, wave, freq, element, level, ST, sep='\t')

    #printing out missing strains per ST, per Wave
    for i in wave_list:
        wave_missing_st_num = (not_st_strains_df[not_st_strains_df['Wave'] == i].shape)[0]
        if wave_missing_st_num > 0:
            wave = "Wave " + str(i)
            gene = element.split('_')[0]
            print(gene, wave, wave_missing_st_num, "No_ST", "No_ST", "No_ST", sep='\t')
            total_unclassified = total_unclassified + wave_missing_st_num

    sum_of_all = total_classified + total_unclassified
    percen_clasified = round(total_classified / sum_of_all * 100,2 )
    percen_unclassified = round(total_unclassified / sum_of_all * 100, 2)

    print(st_num,sum_of_all, total_classified, percen_clasified, total_unclassified, percen_unclassified,sep='\t')
    print()

def main():
    args = parseargs()

    #variables
    start_col = 272
    end_col = 273
    min_strains_per_st = 10
    percen_contam_strains = 15
    min_st_for_figure = 2
    precen_inconsistent_strains = 10
    relate_results_to = "Wave"
    infile_path = '/Users/liamcheneyy/Desktop/vcseventh_22/grapetree/seventh/MGT_isolate_data.txt'

    #df read
    df, want_attributes = df_setup(infile_path, start_col,end_col)

    ##will return MGT level and best ST to describe group
    final_save = {}
    return_df = pd.DataFrame()

    ##calculting non-overlapping STs
    for element in want_attributes:
        print(element)

        ##saving dict
        final_save[element] = {}

        ##create sub_df of element and STs
        col_list = list(df.columns)[0:8] + [element]
        sub_df = df[col_list]

        #calculate the STs for each attribute
        return_dict, return_df, not_st_strains_df = sts_per_attributes(element, sub_df, min_strains_per_st, percen_contam_strains,precen_inconsistent_strains,relate_results_to, df)

        st_freq_and_waves(return_df, return_dict, element, not_st_strains_df, min_st_for_figure)

if __name__ == '__main__':
    main()



#Logic
#to find the STs which best describe attributes

#1. get a list of all strains TRUE and FALSE for gene
#2. then make dict, for each MGT level, get the frequency of STs, the for each ST count the number of TRUE and FALSE for gene
#3. quality filter: remove small STs, remove ST which has too many false and false strains, need just TRUE within threshold
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