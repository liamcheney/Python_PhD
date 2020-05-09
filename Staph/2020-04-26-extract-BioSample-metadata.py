import argparse
from time import sleep as sl
import glob

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args


def main():
    args = parseargs()

    #for extracting metadata from biosamples
    results_dict = {}
    inpath = "/Users/liamcheneyy/Desktop/splits/*xml"
    for file in glob.iglob(inpath):
        infile  = open(file).read().splitlines()

        biosmaple = ""
        collcetion_date = ""
        geo = ""
        flag_country_and_geo = False

        for line in range(0, len(infile), 1):

            if 'namespace="BioSample' in infile[line]:
                biosmaple = infile[line].split('BioSample">')[-1].split('<')[0]


            # if ('collection' in infile[line]) and ('date' in infile[line]):
            #     collcetion_date = infile[line + 1].split('>')[1].strip().split('<')[0]
            #
            if 'geo' in infile[line]:
                geo = infile[line + 1].split('>')[1].strip().split('<')[0]
                flag_country_and_geo = True

            if ('country' in infile[line]) and ("geo" not in infile[line]) and flag_country_and_geo == True:
                geo = infile[line+1].split('>')[1].split('<')[0]

        results_dict[biosmaple] = [geo]

        print(biosmaple,geo)

    with open('/Users/liamcheneyy/Desktop/fixing_saues_metadata.tsv','w') as out:
        for k,v in results_dict.items():
            out.write(k + '\t')
            for i in v:
                out.write(i + '\t')

            out.write('\n')

    #for fixing year metadata
    # infile = open('/Users/liamcheneyy/Desktop/saues_metadata.tsv').read().splitlines()
    #
    # save_dict = {}
    # for line in infile:
    #     col = line.split('\t')
    #     bio = col[0]
    #     date = col[1]
    #
    #
    #     if len(date.split('/')[-1]) == 2:
    #         ##
    #         date = date.split('/')[-1]
    #
    #         if int(date) <= 20:
    #             date = "20" + date
    #             print(bio, date)


    # import pandas as pd
    # ##check bioproject number with metadata
    # infile = open('/Users/liamcheneyy/Desktop/bioprjects_list.csv').read().splitlines()
    #
    # df = pd.read_csv('/Users/liamcheneyy/Desktop/meta_7-gene-alleles.txt', sep='\t')
    #
    # for el in infile:
    #     subdf = df[df['BioProject']==el]
    #     if subdf.shape[0] > 0:
    #         s_subdf = subdf.dropna()
    #         print(el, subdf.shape[0], s_subdf.shape[0])

if __name__ == '__main__':
    main()
