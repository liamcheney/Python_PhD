import argparse
from time import sleep as sl
import os,sys
import time
import multiprocessing as mp

# def parseargs():
#     parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#     parser.add_argument("-i", "--input",
#                          help="Input for annotated contigs.", default="/Users/liamcheneyy/Desktop/out.txt")
#     parser.add_argument("-f", "--fullnamelineage",
#                         help="Input for NCBI Taxonomy fullnamelineage file (/new_taxdump/fullnamelineage.dmp).", default="/Users/liamcheneyy/Downloads/new_taxdump/fullnamelineage.dmp")
#     parser.add_argument("-os", "--output_success",
#                         help="Output for work.",
#                         default="/Users/liamcheneyy/Downloads/alice_GI_worked.txt")
#     parser.add_argument("-of", "--output_failed",
#                         help="Output for work.",
#                         default="/Users/liamcheneyy/Downloads/alice_GI_failed.txt")
#
#     args = parser.parse_args()
#
#     return args
#
# def get_lineage_input(args):
#
#     infile = open(args.input,'r').read().splitlines()
#
#     save_list = []
#     failed_list = []
#     for line in infile:
#         if '[' in line and ']' in line:
#             GI = line.split('\t')[2]
#             lineage = line.split('[')[-1].split(']')[0]
#             save_list.append([GI, lineage])
#         else:
#             failed_list.append(line)
#
#     print("Number of GI read in: " +str(len(save_list)))
#
#     return save_list
#
# def finder(want_lineages, fullname_file, args):
#     for element in want_lineages:
#         lineage = element[-1]
#         GI = element[0]
#         save_list = []
#         for line in fullname_file:
#             if lineage in line:
#                 fix_line = line.split('\t')
#                 fix_line = ';'.join(fix_line).replace(';|','')
#                 save_list.append([GI, lineage, fix_line])
#
#         #multiple hits
#         if len(save_list) > 1:
#             print(line ,len(save_list))
#
# def file_checker(args):
#
#     if os.path.isfile(args.output_success):
#         os.remove(args.output_success)
#
#     with open(args.output_success,'w') as out:
#         pass
#
#
# def main():
#     args = parseargs()
#
#     want_lineages = get_lineage_input(args)
#
#     print("Reading in full lineages.")
#     infile = open(args.fullnamelineage).read().splitlines()
#
#     print('Finding fullnames.')
#     finder(want_lineages, infile, args)

if __name__ == '__main__':
    main()

##previous attempt, way to slow
def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input",
                         help="Input for annotated contigs.", default="/Users/liamcheneyy/Desktop/out.txt")
    parser.add_argument("-f", "--fullnamelineage",
                        help="Input for NCBI Taxonomy fullnamelineage file (/new_taxdump/fullnamelineage.dmp).", default="/Users/liamcheneyy/Downloads/new_taxdump/fullnamelineage.dmp")
    parser.add_argument("-a", "--acc2taxid",
                        help="Input for NCBI Taxonomy accessions2taxid file (/accessions2taxid.zip).",
                        default="/Users/liamcheneyy/Downloads/prot.accession2taxid")
    parser.add_argument("-os", "--output_success",
                        help="Output for work.",
                        default="/Users/liamcheneyy/Downloads/alice_GI_worked.txt")
    parser.add_argument("-of", "--output_failed",
                        help="Output for work.",
                        default="/Users/liamcheneyy/Downloads/alice_GI_failed.txt")

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args

def get_GI_input(args):

    infile = open(args.input,'r').read().splitlines()

    save_list = []
    for line in infile:
        col = line.split('\t')
        GI = col[2]
        save_list.append(GI)

    print("Number of GI read in: " +str(len(save_list)))
    return save_list

def GI_to_taxid(want_GI_list, args):

    ##read in namelienage file
    fullnamelineage_Dict = fullnamelineage_toDict(args)

    #iterate over lagre taxonomy file
    save_dict = mp(want_GI_list, fullnamelineage_Dict, args)

    #check input GI that failed
    failed(want_GI_list, save_dict, args)

def fullnamelineage_toDict(args):

    print('Creating dict for fullenamelineages.')
    ##read in namelienage file
    fullnamelineage = open(args.fullnamelineage, 'r').read().splitlines()

    save_dict= {}
    for line in fullnamelineage:
        col = line.split('	|	')
        taxId = col[0]
        fix_line = line.replace('|',';')
        fix_line = fix_line.replace('\t','')

        save_dict[taxId] = fix_line

    return save_dict

def mp(want_GI_list, fullnamelineage_Dict, args):
    # init objects

    pool = mp.Pool(mp.cpu_count())
    jobs = []

    # create jobs
    print('Read and iterating over large file: accessions2taxid')

    save_dict = {}
    line_count = 0
    with open(args.acc2taxid) as f:
        for line in f:
            jobs.append(pool.apply_async(iterator(line, save_dict, line_count, want_GI_list, fullnamelineage_Dict, args), (line)) )

    # wait for all jobs to finish
    for job in jobs:
        job.get()

    # clean up
    pool.close()

def iterator(line, save_dict, line_count, want_GI_list, fullnamelineage_Dict, args):
    total_lines = 844546480

    for line in f:
            if 'accession' not in line:
                line = line.split('\n')[0]
                col = line.split('\t')
                GI = col[1]
                tax_id = col[2]

                if GI in want_GI_list:
                    lineage = fullnamelineage_Dict[tax_id]
                    save_dict[GI] = lineage
                    print(GI, lineage)
                    with open(args.output_success, 'a') as out:
                            out.write(GI + '\t' + tax_id + '\t' + lineage + '\n')

            pecen_done = (line_count / total_lines ) * 100
            # print(line_count)
            line_count += 1

    percen_found = (len(save_dict.keys()) / len(want_GI_list)) * 100
    print('Found lineages for tax_ids: ' + str(len(save_dict.keys())) + ',' + str(percen_found) + '%')

    write_sucessful(save_dict, args)

def failed(want_GI_list, save_dict, args):
    # check input GI that failed
    failed_list = []
    for el in want_GI_list:
        if el not in save_dict.keys():
            failed_list.append(el)

    ##print failed GIs
    print("Failed input GIs: ")
    print(*failed_list, sep='\t')

    with open(args.output_failed,'w') as out:
        out.write('Faled GIs' + '\n')
        for el in failed_list:
            out.write(el + '\n')

def write_sucessful(save_dict, args):

    with open(args.output_success,'w') as out:
        for key, value in save_dict.items():
            out.write(key + '\t' + value + '\n')

def file_checker(args):

    if os.path.isfile(args.output_success):
        os.remove(args.output_success)

    with open(args.output_success,'w') as out:
        pass

def main():
    args = parseargs()

    file_checker(args)

    ##create list of input GI
    want_GI_list = get_GI_input(args)

    ##iterate over large GI to taxid
    GI_to_taxid(want_GI_list, args)

if __name__ == '__main__':
    main()
