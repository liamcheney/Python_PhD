import argparse
from time import sleep as sl
import pandas as pd

def calculate_per_scheme(df, want_list):

    #make list of MGT ST schemes
    col_list = []
    for column in df:
        # if (str(column.split()[-1]) == 'ST'):
        if column in want_list:
            col_list.append(column)

    #cet frequency of each ST per level
    freq_dict = {}
    avg_dict = {}
    for element in col_list:

        # if element == "MGT4":
        #     sub_df = df[element]
        #     freqs = sub_df.value_counts()
        #     save_freqs = (freqs[freqs > 3])
        #     # print(save_freqs)
        #     save_freqs = save_freqs[save_freqs != 1370]
        #     freq_dict[element] = save_freqs
        #     avg = int(save_freqs.mean())
        #     avg_dict[element] = avg
        #     print(element, avg)
        #
        #     # sub_df.drop(df[element == 9])
        #     # print(sub_df.shape)
        #
        # else:
        sub_df = df[element]
        freqs = sub_df.value_counts()
        save_freqs = (freqs[freqs > 3])
        freq_dict[element] = save_freqs
        avg = int(save_freqs.mean())
        avg_dict[element] = avg
        print(element, avg)

    # for key in
    # calculate average size of

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args


def main():
    args = parseargs()

    want_list = ["MGT2", "MGT3", "MGT4", "MGT5", "MGT6", "MGT7","MGT8", "MGT9"]

    # read_in_metadate_pth = '/Users/liamcheneyy/Downloads/MGT_isolate_data_6.txt'
    read_in_metadate_pth = "/Users/liamcheneyy/Desktop/MGT_isolate_data.txt"

    df = pd.read_csv(read_in_metadate_pth, sep='\t')

    #per lvl dict with size of STs
    st_size_dict = calculate_per_scheme(df, want_list)

if __name__ == '__main__':
    main()
