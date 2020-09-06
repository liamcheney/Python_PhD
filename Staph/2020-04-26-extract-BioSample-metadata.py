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
    abrevs= ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
    args = parseargs()
    #

    #for extracting metadata from a single XML file with multiple biosamples in

    infile = open('/Users/liamcheneyy/Desktop/biosample_result.xml').read()
    files = infile.split('<BioSample access')
    for file in files:
        date = ''
        accession = ''
        age = ''
        sex = ''
        geo = ''
        lines = file.split('   ')
        for line in lines:
            if 'collection date' in line:
                date = line.split('>')[-2].split('<')[0]

            if 'accession' in line:
                accession = line.split('accession')[-1].split('"')[1]

            if 'host age' in line:
                age = line.split('host age')[-1].split('"')[1].strip('>').strip('</Attribute')

            if 'host sex' in line:
                sex = line.split('host sex')[-1].split('"')[1].strip('>').split('<')[0]

            if 'geographic location' in line:
                geo = line.split('geographic location')[-1].split('"')[1].strip('>').split('<')[0]


        print(accession, date, age, sex, geo)

    # #for extracting metadata from biosamples
    # results_dict = {}
    # inpath = "/Users/liamcheneyy/Desktop/splits/*xml"
    # for file in glob.iglob(inpath):
    #     infile = open(file).read().splitlines()
    #     biosmaple = file.split('/')[-1].split('.')[0]
    #     geo = ""
    #     collection_date = ""
    #
    #     for line in range(0, len(infile), 1):
    #
    #         if ('collection' in infile[line]) and ('date' in infile[line]):
    #             collection_date = infile[line + 1].split('>')[1].strip().split('<')[0]
    #
    #             if '-' in collection_date:
    #                 split = collection_date.split('-')
    #                 for x in split:
    #                     if len(x) == 4:
    #                         collection_date = x
    #
    #
    #
    #         # if 'location' in infile[line] or 'country' in infile[line]:
    #         #     geo = infile[line + 1].split('>')[1].strip().split('<')[0]
    #
    #             # print(biosmaple, collection_date)
    #     results_dict[biosmaple] = [collection_date]
    #
    #
    # with open('/Users/liamcheneyy/Desktop/fixing_saues_metadata.tsv','w') as out:
    #     for k,v in results_dict.items():
    #         out.write(k + '\t')
    #         for i in v:
    #             out.write(i + '\t')
    #
    #         out.write('\n')

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


    # ##check bioproject number with metadata
    # import pandas as pd
    # infile = open('/Users/liamcheneyy/Desktop/bioprjects_list.csv').read().splitlines()
    #
    # df = pd.read_csv('/Users/liamcheneyy/Desktop/meta_7-gene-alleles.txt', sep='\t')
    #
    # for el in infile:
    #     subdf = df[df['BioProject']==el]
    #     if subdf.shape[0] > 0:
    #         s_subdf = subdf.dropna()
    #         print(el, subdf.shape[0], s_subdf.shape[0])


    # #taking coverage from assembly metadata
    # for filename in glob.iglob('/Users/liamcheneyy/Downloads/genome_assemblies_asm_stats (2)/ncbi-genomes-2020-05-13/GC*txt'):
    #     infile = open(filename).read().splitlines()
    #     cov = ''
    #     tech = ''
    #     acc = "GCA_" + filename.split('/')[-1].split('_')[1]
    #     date=''
    #
    #     for line in infile:
    #
    #         if "Sequencing technology" in line:
    #             tech = line.split('technology:')[-1].strip('\n')
    #         if "coverage" in line:
    #             cov = line.split(':')[-1].strip('\n')
    #         if "collection" in line:
    #             # cov = line.split(':')[-1].strip('\n')\
    #             print(line)
    #
    #
    #     # print(acc, cov, tech, sep=',')



if __name__ == '__main__':
    main()
