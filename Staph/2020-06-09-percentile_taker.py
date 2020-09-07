import argparse
from time import sleep as sl
import pandas as pd
import statistics
import sys

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args

def get_middle_score(df):
    length=df.shape[0]
    middle = int(length/2)
    scores = sorted(df['Attribute'].values.tolist())
    score = scores[middle]
    print(f'Found middle score: {score}, index: {middle}')
    return middle

def percentile_values(element, l_range, u_range, df):
    want = df.iloc[l_range:u_range, :]
    lowest_val = want['Attribute'].tolist()[0]
    highest_val = want['Attribute'].tolist()[-1]
    print(f'{element}% {lowest_val} {highest_val}', sep='\t')

def get_percentiles(percentiles,middle_index, df):
    save_dict = {}
    for element in percentiles:
        # print(f'Calculating values within {element} percentile of mean.')
        t_length = df.shape[0]
        element_percen = round(float((element / 100) * t_length),2)
        one_tail_percen = int(element_percen / 2)


        l_range = middle_index - one_tail_percen
        u_range = middle_index + one_tail_percen

        ##checking the limits are okay
        if (l_range < 0) or (u_range > t_length):
            if l_range < 0:
                print(f"Only lower limit met: {l_range}")
                sys.exit()
            elif u_range > t_length:
                print(f"Only upper limit met: {u_range}")
                sys.exit()
            else:
                print(f"Both limits met (lower, upper): {l_range} {u_range}")
                sys.exit()

        elif (l_range >= 0) and (u_range <= t_length):
            want = df.iloc[l_range:u_range, :]
            want = want['Locus'].tolist()
            save_dict[element] = want
            percentile_values(element, l_range, u_range, df)

        else:
            print('error')
            sys.exit()

    return save_dict

def add_to_dataframe(percentile_dicts, df):

    for key, value in percentile_dicts.items():
        df[str(key) + "_percentile"] = df['Locus'].isin(value)

    return df

def r_calc_cuttoffs(df):
    save_dict = {}
    infile = open('/Users/liamcheneyy/Desktop/Book7.txt').read().splitlines()
    for line in infile:
        col = line.split('\t')
        peren = col[0]
        lower = float(col[1])
        upper = float(col[2])

        want = df[(df['Attribute'] >= lower) &  (df['Attribute'] <= upper)]
        want = want['Locus'].tolist()
        save_dict[peren] = want

    return save_dict

def main():
    args = parseargs()

    ##will take the percentile from the mean of a dataset
    ##assumes normal distribution for data

    percentiles = [10,20,30,40,50,60,70,80,90,95,96,97,98,99]

    df = pd.read_csv('/Users/liamcheneyy/Desktop/dNdS.txt', sep='\t', dtype={'Attribute':float})
    df = df.sort_values(by=['Attribute'])

    average = round(statistics.mean(df['Attribute']),5)
    median = round(statistics.median(df['Attribute']),5)
    print(f'IMPORTANT: Assumes Normal Distribution')
    print(f"Using Mean: {average}, Median: {median}")
    sl(1)

    #get the middle score
    middle_index = get_middle_score(df)

    #get the loci for percentile at each level
    # percentile_dicts = get_percentiles(percentiles,middle_index, df)

    ##separate function when already have the cut offs
    percentile_dicts = r_calc_cuttoffs(df)

    ##print out
    out_df = add_to_dataframe(percentile_dicts, df)

    out_df.to_csv('/Users/liamcheneyy/Desktop/pref_done.csv', index=False)



if __name__ == '__main__':
    main()
