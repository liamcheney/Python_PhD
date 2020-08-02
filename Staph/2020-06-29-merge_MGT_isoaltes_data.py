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

    MGT_st_df = pd.read_csv('/Users/liamcheneyy/Desktop/MGT_isolate_data (3).txt', sep='\t')
    MGT_st_df = MGT_st_df.drop(columns=['Month','Year','Date','Disease','Host','Type','Source','Postcode','State','Continent','Country'])

    metadata_df = pd.read_excel(open('/Users/liamcheneyy/Desktop/meta_7-gene-alleles.xlsx', 'rb'), sheet_name='Metadata').fillna("none")
    meta_sub = metadata_df[['ID','Year','Country','Continent','BioProject','7 gene ST']]

    merged = MGT_st_df.merge(meta_sub,how='left',left_on='Strain', right_on='ID')

    merged.to_csv('/Users/liamcheneyy/Desktop/meta_MGT_isolate_data (3).txt',sep='\t', index=False)

if __name__ == '__main__':
    main()
