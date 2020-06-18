# import argparse
# from time import sleep as sl
# import os,sys
# import time
# import multiprocessing as mp
#
# ##multi processing
# def process_wrapper(lineID):
#     with open("/Users/liamcheney/Desktop/out.txt") as f:
#         for i,line in enumerate(f):
#             if i != lineID:
#                 continue
#             else:
#                 process(line)
#                 break
#
# #iterating
# def conveter(want_GI_list, fullnamelineage_Dict, args):
#     #init objects
#     core = os.cpu_count()
#     pool = mp.Pool(core)
#     jobs = []
#
#     # create jobs
#     print('Read and iterating over large file: accessions2taxid')
#
#     save_dict = {}
#     line_count = 0
#
#     #create jobs
#     with open(args.acc2taxid) as f:
#         for ID,line in enumerate(f):
#             jobs.append(pool.apply_async(iterator(line, save_dict, line_count, want_GI_list, fullnamelineage_Dict, args), args=(line,)))
#
#     #wait for all jobs to finish
#     for job in jobs:
#         job.get()
#
#     #clean up
#     pool.close()
#
# def iterator(line, save_dict, line_count, want_GI_list, fullnamelineage_Dict, args):
#     total_lines = 844546480
#
#     if 'accession' not in line:
#         line = line.split('\n')[0]
#         col = line.split('\t')
#         GI = col[1]
#         tax_id = col[2]
#
#         if GI in want_GI_list:
#             lineage = fullnamelineage_Dict[tax_id]
#             save_dict[GI] = lineage
#             print(GI, lineage)
#             with open(args.output_success, 'a') as out:
#                     out.write(GI + '\t' + tax_id + '\t' + lineage + '\n')
#
#         pecen_done = (line_count / total_lines ) * 100
#         line_count += 1
#
#     write_sucessful(save_dict, args)
#
# ##storing inforation
# def fullnamelineage_toDict(args):
#
#     print('Creating dict for fullenamelineages.')
#     ##read in namelienage file
#     fullnamelineage = open(args.fullnamelineage, 'r').read().splitlines()
#
#     save_dict= {}
#     for line in fullnamelineage:
#         col = line.split('	|	')
#         taxId = col[0]
#         fix_line = line.replace('|',';')
#         fix_line = fix_line.replace('\t','')
#
#         save_dict[taxId] = fix_line
#
#     return save_dict
# def get_GI_input(args):
#
#     infile = open(args.input,'r').read().splitlines()
#
#     save_list = []
#     for line in infile:
#         col = line.split('\t')
#         GI = col[2]
#         save_list.append(GI)
#
#     print("Number of GI read in: " +str(len(save_list)))
#     return save_list
#
# ##input and output
# def parseargs():
#     parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#     parser.add_argument("-i", "--input",
#                          help="Input for annotated contigs.", default="/Users/liamcheney/Desktop/out.txt")
#     parser.add_argument("-f", "--fullnamelineage",
#                         help="Input for NCBI Taxonomy fullnamelineage file (/new_taxdump/fullnamelineage.dmp).", default="/Users/liamcheney/Downloads/new_taxdump/fullnamelineage.dmp")
#     parser.add_argument("-a", "--acc2taxid",
#                         help="Input for NCBI Taxonomy accessions2taxid file (/accessions2taxid.zip).",
#                         default="/Users/liamcheney/Downloads/sample.txt")
#     parser.add_argument("-os", "--output_success",
#                         help="Output for work.",
#                         default="/Users/liamcheney/Downloads/alice_GI_worked.txt")
#     parser.add_argument("-of", "--output_failed",
#                         help="Output for work.",
#                         default="/Users/liamcheney/Downloads/alice_GI_failed.txt")
#
#     # parser.add_argument("-d", "--database_name", required=True,
#     #                     help="sql database to search (eg. vcseventh)")
#
#     args = parser.parse_args()
#
#     return args
# def file_checker(args):
#
#     if os.path.isfile(args.output_success):
#         os.remove(args.output_success)
#
#     with open(args.output_success,'w') as out:
#         pass
# def failed(want_GI_list, save_dict, args):
#     # check input GI that failed
#     failed_list = []
#     for el in want_GI_list:
#         if el not in save_dict.keys():
#             failed_list.append(el)
#
#     ##print failed GIs
#     print("Failed input GIs: ")
#     print(*failed_list, sep='\t')
#
#     with open(args.output_failed,'w') as out:
#         out.write('Faled GIs' + '\n')
#         for el in failed_list:
#             out.write(el + '\n')
# def write_sucessful(save_dict, args):
#
#     with open(args.output_success,'w') as out:
#         for key, value in save_dict.items():
#             out.write(key + '\t' + value + '\n')
#
# def GI_to_taxid(want_GI_list, args):
#
#     ##read in namelienage file
#     fullnamelineage_Dict = fullnamelineage_toDict(args)
#
#     #iterate over lagre taxonomy file
#     save_dict = conveter(want_GI_list, fullnamelineage_Dict, args)
#
#     # check input GI that failed
#     # failed(want_GI_list, save_dict, args)
#
# def main():
#     args = parseargs()
#
#     file_checker(args)
#
#     ##create list of input GI
#     want_GI_list = get_GI_input(args)
#
#     ##iterate over large GI to taxid
#     GI_to_taxid(want_GI_list, args)
#
# if __name__ == '__main__':
#     main()

##
##trying SQL

import sqlite3
import pandas as pd

##creatinge databasee
mydb = sqlite3.connect("mydatabase")
cursor = mydb.cursor()

##creating a table
# cursor.execute('CREATE TABLE Sample (col1 INT, col2 FLOAT)')

##importing data into database
infile = pd.read_csv('/Users/liamcheney/Desktop/sample.txt',sep='\t', header=0)
# infile.to_sql('Sample', mydb, if_exists="replace")

sub = infile.drop_duplicates()
sub.to_csv('/Users/liamcheney/Desktop/sssample.csv')
##deleting a table
# cursor.execute("DROP table if exists Sample")

##get table names
# cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")

