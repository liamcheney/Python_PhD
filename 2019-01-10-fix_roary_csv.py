##to-date have been writing many individual codes to deal with inconsistencies in strings of nodes between ERR and GCC files
##script will go over roary output and clean up, to avoid later python efforts

from time import sleep as sl
import pandas as pd

infile = '/Users/liamcheneyy/Google Drive/Honours/1_core_genomes/4All_Core_Genome/D2C_i96_99_no_ref_blanks_percen_para_ana.csv'
outfile = '/Users/liamcheneyy/Google Drive/Honours/1_core_genomes/4All_Core_Genome/D2C_i96_99_edited_names_percen_para_ana.csv'
df = pd.read_csv(infile, sep=',')

##fixing length
for column_head in df:
    if 'ERR' in column_head:
        df[column_head] = df[column_head].str.replace('length_','')

#fixing the GC columns
col_head_list = []
for column_head in df:
    if 'GC' not in column_head:
        col_head_list.append(column_head)
    if 'GC' in column_head:
        new_header = 'GCA' + str(column_head[4:])
        col_head_list.append(new_header)
        df[column_head] = df[column_head].str.replace('GCA_','GCA')
df.columns = col_head_list

#fixing little stuffs
for column_head in df:
    if 'ERR' in column_head:
        df[column_head] = df[column_head].str.replace('NODE_', 'NODE')
    if 'GC' in column_head:
        df[column_head] = df[column_head].str.replace('NZ_', 'NZ')
    if 'GC' in column_head:
        df[column_head] = df[column_head].str.replace('NC_', 'NC')

df.to_csv(outfile, sep=',', index=False)


# #fixing the node names etc in the FNA files
# from Bio import SeqIO
# import glob
#
# infile_path = '/Users/liamcheneyy/Desktop/untitled folder/'
#
# for filename in glob.iglob(infile_path + 'ERR*'):
#     file_out = filename.split('/')[-1]
#     outfile = open('/Users/liamcheneyy/Desktop/untitledfolder2/' + file_out, 'w')
#     ###fixing ERR
#     for record in SeqIO.parse(filename, 'fasta'):
#         record.id = str(record.id).split('_')[1]
#         outfile.write(str('>' + record.id + '\n'))
#         outfile.write(str(record.seq + '\n'))
    ###fixing GCC
    # for record in SeqIO.parse(filename, 'fasta'):
    #     if len(record.id.split('_')) > 1:
    #         record.id = str(record.id).split('_')[0] + str(record.id).split('_')[1]
    #         outfile.write(str('>' + record.id + '\n'))
    #         outfile.write(str(record.seq + '\n'))
    #     else:
    #         print(filename)
    #         outfile.write(str('>' + record.id + '\n'))
    #         outfile.write(str(record.seq + '\n'))
    # outfile.close()

