import argparse
from time import sleep as sl
import pandas as pd

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args


def main():
    args = parseargs()


    infile_path = "/Users/liamcheneyy/Desktop/continent_specific_STs.csv"
    df = pd.read_csv(infile_path)

    contamination_amount = 20
    min_size_st = 40

    #get list of STs
    MGT3_STs = df['MGT3'].unique().tolist()

    save_dict = {}
    #check the continents for each MGT3 ST

    for i in MGT3_STs:
        sub_df = df[df['MGT3'] == i]
        if list(sub_df.shape)[0] > min_size_st:
            continents_count_df = sub_df['Continent'].value_counts().to_dict()
            contam_total = sum(sub_df['Continent'].value_counts().tolist()[1:])
            total = sum(sub_df['Continent'].value_counts().tolist())

            first_continent_count = continents_count_df[list(continents_count_df.keys())[0]]
            first_continent = list(continents_count_df.keys())[0]
            if (contam_total / total * 100) <= contamination_amount:
                print("CONT SPECIFIC","ST" + str(i), first_continent_count, first_continent, continents_count_df,  sep='\t')
            else:
                print("MULTI CONT","ST" + str(i), first_continent_count, "", continents_count_df,  sep='\t')



if __name__ == '__main__':
    main()
