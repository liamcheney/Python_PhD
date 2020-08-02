import argparse
from time import sleep as sl
import multiprocessing as mp
import pandas as pd

def strains_w_both_meta(df):

    sub = df[['Country','Year']].dropna()
    new_count = sub.shape[0]

    return new_count

def min_max_year(df):
    sub = df[['Country','Year']].dropna()
    sub = sub.sort_values('Year', ascending=True)

    sub_year = list(sub['Year'])
    min = int(sub_year[0])
    max = int(sub_year[-1])

    return min, max

def unique_list(column):
    save_list = []
    save_list = list(column.unique())

    return save_list

def get_STs_w_metadata(df):

    unique_ST_list = unique_list(df['7 gene ST'])

    for ST in unique_ST_list:
        sub = df[df['7 gene ST'] == ST]

        ST_size = sub.shape[0]
        if ST_size > 50:

            ST_genome_w_both_meta = strains_w_both_meta(sub)

            ST_genome_w_year = sub['Year'].dropna().shape[0]
            ST_genome_w_country = sub['Country'].dropna().shape[0]

            percen_w_meta = round(ST_genome_w_both_meta / ST_size *100, 1)
            print(ST, ST_size, ST_genome_w_both_meta, percen_w_meta, ST_genome_w_year, ST_genome_w_country, sep='\t')

def get_Country_w_metadata(df):

    unique_country_list = unique_list(df['Country'])
    for country in unique_country_list[1:]:
        sub = df[df['Country']==country]

        # #how many genomes per country with time
        country_genome_with_time_meta = strains_w_both_meta(sub)
        country_size = sub.shape[0]
        percen_total = round(country_genome_with_time_meta / country_size * 100, 1)

        #number of STs for a country
        unique_STs = unique_list(sub['7 gene ST'])
        num_unique_STs = len(unique_STs)

        if country_genome_with_time_meta > 0:
            min_year, max_year = min_max_year(sub)
            year_range = str(min_year) + '-'+ str(max_year)
        else:
            year_range = 0


        print(country, country_size, country_genome_with_time_meta,percen_total, num_unique_STs, year_range, sep='\t')

if __name__ == '__main__':

    df = pd.read_csv('/Users/liamcheney/Desktop/test.txt',sep='\t')

    # both_meta_count = strains_w_both_meta(df)

    # STs_w_metadata = get_STs_w_metadata(df)
    Country_w_metadata = get_Country_w_metadata(df)
