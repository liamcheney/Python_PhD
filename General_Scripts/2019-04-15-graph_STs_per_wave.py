from time import sleep as sl
import pandas as pd

infile_path = "/Users/liamcheneyy/Desktop/vcseventh_15/grapetree/all_metadata.txt"
list_of_levels = ['MGT2 ST', 'MGT3 ST', 'MGT4 ST']
out_path = "/Users/liamcheneyy/Desktop/MGT5_STs_AMRs.tsv"

def read_in(infile_path):

    df = pd.read_csv(infile_path, sep='\t')

    wave1_df = df[df['Wave'] == 1.0]
    wave2_df = df[df['Wave'] == 2.0]
    wave3_df = df[df['Wave'] == 3.0]

    return df, wave1_df, wave2_df, wave3_df

##for waves and ST
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
def wave_to_st(infile_path, list_of_levels, out_path):

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

##for ST and AMRS all waves
def STs_and_AMRs(infile_path):

    MGT_lvl_want = 'MGT5'
    st_number_choose = 0

    df, wave1_df, wave2_df, wave3_df = read_in(infile_path)

    list_of_MGT4_sts = df[MGT_lvl_want + ' ST'].value_counts()
    list_of_MGT4_sts = list_of_MGT4_sts[list_of_MGT4_sts > st_number_choose]

    column_headers = ["ST" + str(x) for x,i in list_of_MGT4_sts.iteritems()]
    print(column_headers)
    result_df = pd.DataFrame(columns=column_headers)
    for label,value in list_of_MGT4_sts.iteritems():
        temp_df = df[df[MGT_lvl_want + ' ST'] == label]
        df_AMR_count = temp_df.iloc[:,59:].apply(lambda x: x.count(), axis=0)
        pecen_AMR_count = df_AMR_count.apply(lambda x: float((x/value)*100))
        result_df['ST' + str(label)] = pecen_AMR_count

    result_df = result_df.round(2)
    
    result_df.to_csv('/Users/liamcheneyy/Desktop/' + MGT_lvl_want + '_ST_AMRs_a.csv',sep=',')

##for waves and AMRS
def wave_to_st(infile_path, list_of_levels, out_path):

    df, wave1_df, wave2_df, wave3_df = read_in(infile_path)

    wave1_MAR_count = wave1_df.iloc[:,59:].apply(lambda x: x.count(), axis=0)
    # wave1_MAR_count = wave1_MAR_count.sort_values(ascending=False)

    wave2_MAR_count = wave2_df.iloc[:,59:].apply(lambda x: x.count(), axis=0)
    # wave2_MAR_count = wave2_MAR_count.sort_values(ascending=False)


    wave3_MAR_count = wave3_df.iloc[:,59:].apply(lambda x: x.count(), axis=0)
    # wave3_MAR_count = wave3_MAR_count.sort_values(ascending=False)

    wave1_MAR_count.to_csv('/Users/liamcheneyy/Desktop//wave1_MAR_count.csv',sep=',')
    wave2_MAR_count.to_csv('/Users/liamcheneyy/Desktop//wave2_MAR_count.csv',sep=',')
    wave3_MAR_count.to_csv('/Users/liamcheneyy/Desktop//wave3_MAR_count.csv',sep=',')



def main(infile_path, list_of_levels, out_path):

    #if want sts in each wave
    STs_and_AMRs(infile_path)

    #finding antibiotic resistance genes
    # wave_to_st(infile_path, list_of_levels, out_path)

main(infile_path, list_of_levels, out_path)