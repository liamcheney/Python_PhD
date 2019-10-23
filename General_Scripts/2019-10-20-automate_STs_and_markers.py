import argparse
from time import sleep as sl
import pandas as pd

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args

def sts_per_attributes(element, df):

    #create sub dataframe of MGT ST levels and attribute
    sub_df_cols = [element] + [x for x in df if 'MGT' in x and 'CC' not in x and 'MGT9' not in x]
    sub_df = df[sub_df_cols]
    mgt_col_list = [x for x in sub_df if "MGT" in x]

    for item in mgt_col_list:
        want_levels = [sub_df_cols[0], item]
        sub_level_df = sub_df[want_levels]
        st_freq_dict = sub_level_df[item].value_counts()
        for index, row in sub_level_df.iterrows():
            row_values = row.tolist()

            sl(1)




def main():
    args = parseargs()

    #variables
    start_col = 28
    end_col = 29
    include_CCs = False

    #read in metadata
    df = pd.read_csv('/Users/liamcheneyy/Desktop/vcseventh_21/grapetree/seventh/metadata.csv', index_col='ID', low_memory=False)

    #columns of attributes to take
    want_attributes = list(df.columns.values[start_col:end_col])

    #calculate the STs for each attribute
    save_dict = {}
    for element in want_attributes:
        save_dict[element] = sts_per_attributes(element, df)


if __name__ == '__main__':
    main()

#Logic
#to find the STs which best describe attributes

#challenges
#2. the script must choose the best MGT level
    #find the level which defines the least amount of STs, with the largest number.
#3. script must choose an ST unique to a certain attribute
