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

from multiprocessing import Pool

def GI_getter(infile_path):

    infile = open(infile_path).read().splitlines()

    GI_list = []
    beg_list = []
    for line in infile:
        col = line.split('\t')
        GI = col[2]
        beg = GI[0:3]

        if GI not in GI_list:
           GI_list.append(GI)

        if beg not in beg_list:
            beg_list.append(beg)

    return GI_list, beg_list

def gi_to_taxID_dict_creator(GI_list, beg_list, pro_2_tax_path):

    print('Converting GIs to taxIDs. Can take a while.')
    save_dict = {}

    all_GI_len = len(GI_list)

    ##go over the inputs for beg list
    count = 1
    failed_list = []
    for beg in beg_list:

        if os.path.isfile(pro_2_tax_path + '/' + beg):
            with open(pro_2_tax_path + '/' + beg) as f:
                for lines in f:
                    line = lines.splitlines()[0]
                    col = line.split('\t')
                    in_GI = col[1]
                    if in_GI in GI_list:
                        tax_id = col[2]
                        save_dict[in_GI] = {'taxid':tax_id, 'lineage':''}
                        print(f"Found {count} out of {all_GI_len}. Database {beg}. {in_GI}")
                        count += 1
        else:
            failed_list.append(beg)

    return save_dict, failed_list

def taxID_to_lineage_creator(GI_to_taxID_dict, full_name):

    print('Getting full lineage names.')

    ##conver to full lineage
    save_dict = {}
    with open(full_name) as f:
        for lines in f:
            line = lines.splitlines()[0]
            tax_id = line.split('\t')[0]
            lineage = ' '.join(line.replace('|','').replace('\t',' ').split(' ')[1:]).strip(' ')

            for key, value in GI_to_taxID_dict.items():
                if tax_id in value['taxid']:
                    GI_to_taxID_dict[key]['lineage'] = lineage

    return GI_to_taxID_dict

def iterator(GI_list, beg_list, full_name, pro_2_tax_path):

    GI_to_taxID_dict, failed_list = gi_to_taxID_dict_creator(GI_list, beg_list, pro_2_tax_path)

    taxID_to_lineage_dict = taxID_to_lineage_creator(GI_to_taxID_dict, full_name)

    return taxID_to_lineage_dict, failed_list

def outwrite(GI2lineage_dict, infile_path, output_path, failed_list):

    print(f"Writing out results.")
    with open(output_path,'w') as out:
        infile = open(infile_path).read().splitlines()
        for line in infile:
            col = line.split('\t')
            GI = col[2]

            taxid = GI2lineage_dict[GI]['taxid']
            lineage = GI2lineage_dict[GI]['lineage']

            out.write(line + '\t' + taxid + '\t' + lineage + '\n')

    print(f"{len(failed_list)} failed as there were not protein2taxid databases for: {failed_list}")

def main():
    starttime = time.time()

    ##inputs
    infile_path = '/Users/liamcheneyy/Desktop/out.txt'
    full_name = '/Users/liamcheneyy/Desktop/fullnamelineage.dmp'
    pro_2_tax_path = '/Users/liamcheneyy/Desktop/split_protaccessions2taxid/'
    output_path = '/Users/liamcheneyy/Desktop/out_test.txt'

    ##get list of accessions
    GI_list, beg_list = GI_getter(infile_path)

    ##iterator
    GI2lineage_dict, failed_list = iterator(GI_list, beg_list, full_name, pro_2_tax_path)

    ##write out
    outwrite(GI2lineage_dict, infile_path, output_path, failed_list)

    endtime = time.time() - starttime
    print(f"Finished in {endtime}")

if __name__ == "__main__":
    main()