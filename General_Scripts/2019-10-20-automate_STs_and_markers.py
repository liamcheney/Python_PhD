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

    sub_df_cols = [element] + [x for x in df if 'MGT' in x and 'CC' not in x and 'MGT9' not in x]
    sub_df = df[sub_df_cols]
    print(sub_df)



def main():
    args = parseargs()

    #variables
    start_col = 28
    end_col = 29
    include_CCs = False

    #read in metadata
    df = pd.read_csv('/Users/liamcheneyy/Desktop/metadata.csv', index_col='ID', low_memory=False)

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
#3. script must choose an ST unique to a certain attribute
