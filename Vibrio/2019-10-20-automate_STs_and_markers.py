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
    df = pd.read_csv(infile_path, index_col='ID', low_memory=False, sep='\t')

    df[list(df.columns)[0:8]] = df[list(df.columns)[0:8]].astype(int)

    #columns of attributes to take
    want_attributes = list(df.columns.values[start_col:end_col])

    return df, want_attributes
def sts_per_attributes(element, sub_df, min_strains_per_st, percen_contam_strains,precen_inconsistent_strains, df):
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
    STs_list = list(st_freq_dict.keys())
    final_ST_list = []
    for x in STs_list:
        y = x.replace(' ','_')
        final_ST_list.append(y)
    
    print(f'{element} {num_of_st} {all_strains} {percen_all_strains} {total_strains_num} {number_TRUE} {percen_with_described} {non_described_strains} {final_ST_list}')
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

def main():
    args = parseargs()

    #variables
    start_col = 22
    end_col = 32
    min_strains_per_st = 10
    percen_contam_strains = 20
    precen_inconsistent_strains = 10
    infile_path = '/Users/liamcheney/Desktop/newgenotyping/genotyping_w_meta.txt'

    #df read
    df, want_attributes = df_setup(infile_path, start_col,end_col)

    ##calculting non-overlapping STs
    print(f'Genotype #ST Total_Strains %ofStrains #wGenotype #Described %Described #NotDescribed STs')
    for element in want_attributes:
        ##create sub_df of element and STs
        col_list = list(df.columns)[1:8] + [element]
        sub_df = df[col_list]

        #calculate the STs for each attribute
        sts_per_attributes(element, sub_df, min_strains_per_st, percen_contam_strains,precen_inconsistent_strains, df)

if __name__ == '__main__':
    main()