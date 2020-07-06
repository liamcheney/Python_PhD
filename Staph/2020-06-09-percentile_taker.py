import argparse
from time import sleep as sl
import pandas as pd

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args

def get_middle_index(df, average_atr):

    df_sort = df.iloc[(df['Attribute'] - average_atr).abs().argsort()[:2]]
    df_sort = df_sort.index.tolist()
    middle_avgerage_index =df_sort[0] + 1
    return middle_avgerage_index

def get_percentiles(percentiles,middle_avgerage_index, df):
    save_dict = {}
    for element in percentiles:
        t_length = df.shape[0]
        percen_length = float((element / 100) * t_length)
        one_tail_percen_length = float(percen_length / 2)

        l_range = middle_avgerage_index - one_tail_percen_length
        u_range = middle_avgerage_index + one_tail_percen_length
        if (l_range < 0) or (u_range > t_length):
            print("Percentile met limits (lower, upper): ", l_range, u_range, sep='\t')

            if l_range < 0:
                print("lower limit met")

            elif u_range > t_length:
                print("upper limit met")
            l_range = 0
            u_range = t_length

        take_l_range = int(l_range)
        take_u_range = int(u_range)

        want = df.iloc[take_l_range:take_u_range, :]
        want = want['Locus'].tolist()
        # print(len(want))
        save_dict[element] = want

    return save_dict

def add_to_dataframe(percentile_dicts, df):

    for key, value in percentile_dicts.items():
        df[str(key) + "_percentile"] = df['Locus'].isin(value)

    return df

def main():
    args = parseargs()

    ##will take the percentile from the mean of a dataset
    ##assumes normal distribution for data

    percentiles = [95,96]

    df = pd.read_csv('/Users/liamcheney/Desktop/Book2.txt', sep='\t', dtype={'Attribute':float})

    average_atr = float(pd.DataFrame.mean(df['Attribute']))
    middle_avgerage_index = get_middle_index(df, average_atr)
    print(f"Using average of: {average_atr}")

    percentile_dicts = get_percentiles(percentiles,middle_avgerage_index, df)

    add_to_dataframe(percentile_dicts, df)

    df.to_csv('/Users/liamcheney/Desktop/pref_done.csv', index=False)



if __name__ == '__main__':
    main()
