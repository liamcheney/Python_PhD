import argparse
from time import sleep as sl
import glob
from Bio import SeqIO

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args


def main():
    args = parseargs()

    ###COUNTING GENES PER STRAIN
    # infiles = open('/Users/liamcheneyy/Desktop/all_resfinder_outputs.csv').read().split('\n\n')
    # infiles = open('/Users/liamcheneyy/Desktop/all.csv').read().split('./')
    #
    # print(len(infiles))
    # save_dict=  {}
    # for x in infiles[1:]:
    #     infile = x.splitlines()
    #     num = len(infile) - 2
    #     acc = infile[0].strip('_out.csv').split('/')[-1]
    #     save_dict[acc] = num
    #
    # with open('/Users/liamcheneyy/Desktop/all_out.csv','w') as out:
    #     for k, v in save_dict.items():
    #         out.write(k + ',' + str(v) + '\n')


    # ##COUNTING FREQUENCIES OF GENES
    # genes_list = open('/Users/liamcheneyy/Desktop/card_genes.txt').read().splitlines()
    # gene_dict = {}
    # for line in genes_list:
    #     gene_dict[line] = 0
    #
    # infiles = open('/Users/liamcheneyy/Desktop/all_card.csv').read().split('./')
    #
    # for file in infiles[1:-1]:
    #     infile = file.splitlines()
    #     if len(infile) > 2:
    #         for line in infile[2:]:
    #             col = line.split(',')
    #             gene = col[4]
    #             gene_dict[gene] = gene_dict[gene] + 1
    #
    # with open('/Users/liamcheneyy/Desktop/all_out.csv','w') as out:
    #     for k, v in gene_dict.items():
    #         out.write(k + ',' + str(v) + '\n')

    # ##MOST COMMON COMBINATIONS

    combination_dict ={}

    infiles = open('/Users/liamcheneyy/Desktop/all_card.csv').read().split('./')

    for file in infiles[1:-1]:
        infile = file.splitlines()
        acc = infile[0].strip('_out.csv').split('/')[-1]

        save = []
        if len(infile) > 2:
            for line in infile[2:]:
                col = line.split(',')
                gene = col[4]

                ##check for dups
                if gene not in save:
                    save.append(gene)

        save = sorted(save)
        save_str = '/'.join(save)

        if save_str not in combination_dict.keys():
            combination_dict[save_str] = 1

        if save_str in combination_dict.keys():
            combination_dict[save_str] = combination_dict[save_str] + 1

    #check that duplicates with different orders dont exist
    list1 = []
    list2 = []

    for i in combination_dict.keys():
        splited = i.split('/')
        list1.append(splited)
        list2.append(splited)

    for q in list1:
        count = 0
        for l in list2:
            if set(q) == set(l):
                count = count + 1

        if count > 1:
            print(q)


    # with open('/Users/liamcheneyy/Desktop/all_out.csv','w') as out:
    #     for k, v in combination_dict.items():
    #         out.write(k + ',' + str(v) + '\n')


if __name__ == '__main__':
    main()

##wanted
#number of plasmids
