import argparse
from time import sleep as sl
import pandas as pd

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    args = parser.parse_args()

    return args


def main():
    args = parseargs()

    df = pd.read_csv('/Users/liamcheneyy/Desktop/metadata.txt', sep='\t', na_values=' ')
    # df = df[((df['Continent'] == 'AFRICA')) & (df['Wave'] == 3) ]
    df_number_of_st = df.groupby('MGT5 ST').size()
    df_number_of_st = df_number_of_st.sort_values(ascending=False)
    number_of_sts = len(df_number_of_st)



    print(df_number_of_st)

    st_total = 0
    check_st_total = 0
    describe_list = []
    #check if ST is found over threshold of strains
    for row, value in df_number_of_st.iteritems():
        st_total = st_total + int(value)

        if int(value) > 5 :
            check_st_total = check_st_total + int(value)
            describe_list.append(row)
            # print(row,value)

    print(st_total)
    print(check_st_total)
    print(check_st_total / st_total * 100)
    print(sorted(describe_list))


if __name__ == '__main__':
    main()
