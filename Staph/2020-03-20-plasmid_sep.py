import argparse
from time import sleep as sl
from Bio import SeqIO
import glob

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args


def main():
    args = parseargs()

    # infiles = open("/Users/liamcheney/Desktop/all_ncbi_outputs.csv").read().split('\n\n')
    #
    # save_dict = {}
    # save_list = []
    # for file in infiles:
    #     infile = file.split('\n')
    #     for line in infile[1:]:
    #         col = line.split(',')
    #         gene = col[4]
    #         if gene not in save_list:
    #             save_list.append(gene)

    # print(len(save_list))
    # for i in save_list:
    #     print(i)

    # save_dict = {}
    # for filename in glob.iglob('/Users/liamcheney/Desktop/outs/*_out.txt'):
    #     infile = open(filename).read().splitlines()
    #     acc = filename.split('/')[-1].split('_')[0]
    #     intact = 0
    #     partial = 0
    #     called = 0
    #     print(acc)
    #
    #     for line in infile:
    #
    #         if "intact:" in line:
    #             intact = int(line.split(':')[-1].strip())
    #         if "partial:" in line:
    #             partial = int(line.split(':')[-1].strip())
    #         if "called:" in line:
    #             called = int(line.split(':')[-1].strip())
    #
    #     total_pass = intact + partial + called
    #     failed =  1860 - total_pass
    #     save_dict[acc] = {'total':total_pass, 'fail':failed}
    #
    # with open('/Users/liamcheney/Desktop/all_out11.csv','w') as out:
    #     for k,v in save_dict.items():
    #         out.write(k + ',')
    #         for x in v:
    #             out.write(str(save_dict[k][x]) + ',')
    #         out.write('\n')


    ##check how gene went across dataset
    infile = open('//Users/liamcheney/Desktop/all_alleles_out.txt').read().splitlines()

    #make list of genes
    save_dict = {}
    for line in infile[0:1]:
        col = line.split(' ')
       # acc=col[0].split('/')[-1].split('_')[0]
        for j in col[4:]:
            gene = j.split(':')[0].strip('>')
            save_dict[gene] = {'pass':0, 'fail':0}

    #iterate over
    for line in infile:
        acc=col[0].split('/')[-1].split('_')[0]
        col = line.split(' ')
        for j in col[4:]:
            gene = j.split(':')[0].strip('>')
            if ':0' in j:
                save_dict[gene]['fail'] = save_dict[gene]['fail'] + 1
            else:
                save_dict[gene]['pass'] = save_dict[gene]['pass'] + 1

    for k,v in save_dict.items():
        print(k,v['pass'], ['fail'])

if __name__ == '__main__':
    main()