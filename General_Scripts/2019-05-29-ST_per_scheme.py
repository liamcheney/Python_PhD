import argparse
from time import sleep as sl
from pandas import ExcelWriter
import pandas as pd
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt


def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-i", "--input", required=True,
                        help="csv input of STs")
    parser.add_argument("-o", "--output", required=True,
                        help="output file .csv")

    args = parser.parse_args()

    return args


def get_sts_per_scheme(args):

    df = pd.read_csv(args.input, sep=',')
    keep_col_list = [x for x in df.columns.values if 'MGT' in x and 'DST' not in x]
    df = df[df.columns.intersection(keep_col_list)]

    result_dict = {}
    for col in df:
        value_count = df[col].value_counts()
        result_dict[col] = value_count

    xls_path = args.output
    writer = pd.ExcelWriter(xls_path, engine='xlsxwriter')

    for key, value in result_dict.items():
        result_dict[key].to_excel(writer, sheet_name=str(key))
        print(len(value))
    writer.save()


def main():
    args = parseargs()
    result_dict = get_sts_per_scheme(args)

if __name__ == '__main__':
    main()
