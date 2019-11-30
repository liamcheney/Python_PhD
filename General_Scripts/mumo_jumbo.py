# import glob
# from time import sleep as sl
# import matplotlib
# # matplotlib.use("TkAgg")
# # from matplotlib import pyplot as plt
# # from collections import Counter
# #
# # infile_path = '/Users/liamcheneyy/Desktop/failure_reasons/'
#
# # genome_size_list = []
# # large_contig_list = []
# # kraken_list = []
# # n50_list = []
# #
# # contam_count = 0
# # gen_count = 0
# # for filename in glob.iglob(infile_path + '*'):
# #     infile = open(filename).read().splitlines()
# #     file_acc = filename.split('/')[-1].strip('_failure_reason.txt')
# #     for line in infile:
# #         if "Genome outside of allowed range:" in line:
# #             genome_size = line.split('(')[-1].strip(')')
# #             genome_size_list.append([file_acc,float(genome_size)])
# #             # gen_count = gen_count + 1
# #
# #         if "largest contig <" in line:
# #             l_contig = line.split('(')[-1].strip(')')
# #             large_contig_list.append(int(l_contig))
# #             gen_count = gen_count + 1
# #
# #
# #         if "N50 less" in line:
# #             n50 = int(line.split('(')[-1].strip(')')) * 100
# #             n50_list.append(int(n50))
# #             # gen_count = gen_count + 1
# #
# #
# #         if "contaminated" in line:
# #             contam_species = line.split('\t')[-1]
# #             kraken_list.append(contam_species)
# #             contam_count = contam_count + 1
# #
# # for i,b in genome_size_list:
# #
# #     print(i)
# # plt.hist(kraken_list, bins=40)
# # plt.xticks(rotation=90)
# # plt.show()
# #
# # kraken = Counter(kraken_list)
# #
# # import plotly.plotly as py
# # import plotly.graph_objs as go
# #
# # import pandas as pd
# #
# # df = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/2014_us_cities.csv')
# #
# # df['text'] = df['name'] + '<br>Population ' + (df['pop']/1e6).astype(str)+' million'
# # limits = [(0,2),(3,10),(11,20),(21,50),(50,3000)]
# # colors = ["rgb(0,116,217)","rgb(255,65,54)","rgb(133,20,75)","rgb(255,133,27)","lightgrey"]
# # cities = []
# # scale = 5000
# # #
# # # for i in range(len(limits)):
# # #     lim = limits[i]
# # #     df_sub = df[lim[0]:lim[1]]
# # #     city = go.Scattergeo(
# # #         locationmode = 'USA-states',
# # #         lon = df_sub['lon'],
# # #         lat = df_sub['lat'],
# # #         text = df_sub['text'],
# # #         marker = go.scattergeo.Marker(
# # #             size = df_sub['pop']/scale,
# # #             color = colors[i],
# # #             line = go.scattergeo.marker.Line(
# # #                 width=0.5, color='rgb(40,40,40)'
# # #             ),
# # #             sizemode = 'area'
# # #         ),
# # #         name = '{0} - {1}'.format(lim[0],lim[1]) )
# # #     cities.append(city)
# # #
# # # layout = go.Layout(
# # #         title = go.layout.Title(
# # #             text = '2014 US city populations<br>(Click legend to toggle traces)'
# # #         ),
# # #         showlegend = True,
# # #         geo = go.layout.Geo(
# # #             scope = 'usa',
# # #             projection = go.layout.geo.Projection(
# # #                 type='albers usa'
# # #             ),
# # #             showland = True,
# # #             landcolor = 'rgb(217, 217, 217)',
# # #             subunitwidth=1,
# # #             countrywidth=1,
# # #             subunitcolor="rgb(255, 255, 255)",
# # #             countrycolor="rgb(255, 255, 255)"
# # #         )
# # #     )
# # #
# # # fig = go.Figure(data=cities, layout=layout)
# # # py.iplot(fig, filename='d3-bubble-map-populations')
#
# # infile_path = open('/Users/liamcheneyy/Desktop/Untitled.txt').read().splitlines()
# #
# # issue_list = ["virulence","chemotaxis","biotic","pilus",'membrane','secretion','t6ss','t1ss','t2ss','toxin','periplasmic','toxin','periplasmic','extracellular','capsid','capsular','wall'
# #               ,'flagella','mobile','plasmid','pilin','drug','receptor','resistance','peptidoglycan','transposase','twin','fimbrial','transporter']
# #
# # non_issue_dict = {}
# # issue_dict = {}
# # for line in infile_path:
# #     cell = line.split('\t')
# #     for el in issue_list:
# #         if el in cell[1].lower():
# #             issue_dict[cell[0]] = cell[1]
# #         else:
# #             non_issue_dict[cell[0]] = cell[1]
# #
# # # print(len(issue_dict.keys()))
# # for key, value in issue_dict.items():
# #     print(key + '\t' + value)
#
# # infile = open('/Users/liamcheneyy/Desktop/vc_core_genes.txt','r').read().splitlines()
# # bio_list= list([x for x in open('/Users/liamcheneyy/Desktop/biocycle_list.txt').read().splitlines()])
# #
# # for line in infile:
# #     if line in bio_list:
# #         print(line + '\t' + "TRUE")
# #     else:
# #         print(line + '\t' + "FALSE")
#
# # bio_list= list([x for x in open('/Users/liamcheneyy/Desktop/biocycle_list.txt').read().splitlines()])
#
# # for i in open('/Users/liamcheneyy/Desktop/alice_roary_dict_with_cord.txt').read().splitlines():
# #     col = i.split('{')
# #     for j in col:
# #         print(j)
#
# # infile = open('/Users/liamcheneyy/Desktop/Untitled.txt','r').read().splitlines()
# # new_list = []
# # for lin in infile:
# #     if lin != "":
# #         if '\t' in lin:
# #             frag_list = []
# #             for frag in lin.split('\t'):
# #                 col = frag.split('_')
# #                 gene_id = col[1]
# #                 new = "CS_01_" + gene_id
# #                 frag_list.append(new)
# #             final_new = '\t'.join(frag_list)
# #             print(final_new)
# #
# #         elif '\t' not in lin:
# #             col = lin.split('_')
# #             gene_id = col[1]
# #             new = "CS_01_" + gene_id
# #             print(new)
# #     else:
# #         print(lin)
# #
# # # for line in infile:
# # #     col = line.split('\t')
# # #     if col[1] in my_strains:
# # #         if col[2] != '/':
# # #             print(col[1:])
# #
# #
# # i="atgttttcaaaggagtcgttacaaagcttctcgatcactaagaaggatgccctttatgat\
# # cttcttaatcgggaagatctggattggcgcaagttgcaccttcaaacagcgaagaccgtc\
# # atcagtttgactgaccacagcgccatacgcgcctatgttgttgatgattctgttaaaaca\
# # cggcgcggtaaaaccatgccaggcatctctagccactttgatcacctcaatggtcgttgt\
# # gtgatggggcagcaaatcctgaccctagggctggcgactgaaaagcaatttatcccttta\
# # gatagcgaactctatatcagtcgagttaaatcacaacctctgaccaaagcctttgccgat\
# # ggtcgcagtatcgccgcaaaacgttatcgtgatgcgcagtcaatgactaaaccggagatg\
# # gtccatggcatgattaagcgagccgatcggtccggcattcatgcacagtattttttagcc\
# # gattcatggtttgccagtaagtccatgctggcctttatggaagcgcaatcgttagtctcc\
# # attgtccgtatgaagaaaaataaaatgagttaccgcgtggccggtagtgacaaggttctc\
# # tcaacagcggctgagctttaccaacaccacatcaaaggccaatggcagaagattaagggt\
# # agaccataccagtctaaggcgatcaccgttgagttgaatctagcaaaaagcatcaaagag\
# # cctgatcactggatcaaggtcaagctgttgtttgttcgtggtgtcaacgaagaaaaacag\
# # cgcgctggaaagcatgattgggccttgtttttatccaccgacactcacctcagtgatgag\
# # agaatactcgaaatctatgcgttgcgctggggcatcgaagtgtactttaaggaagcgaaa\
# # caaaagcttgggttcctcaaagagcaaagtacacattacagtgcttatatcgcgtccatc\
# # catttaacggcgctgcgattctgcttgctgctactgactcaacacgaggaaggtgccgct\
# # agactgagtgatagtcgtaacgacatgatcaacagtctgtgcaccttagattttgctagc\
# # cgactctgggttatctttagagcactgatatcgggagcgcttgatgagctcagcaagcta\
# # tacggagtcagcgccgcacaagagatcatgaaccaaatagacaaaacggttcaggagttc\
# # tttatgcaagtcatgcaaatggataccttcacgctgaggctagaagccaaatctattggt\
# # gatgagtgctga"
# #
# # print(len(i))
# # with open('/Users/liamcheneyy/Desktop/test.txt','w') as out:
# #     out.write(str(i.upper()))
#
# from Bio import SeqIO
# import glob
#
# # infile = '/Users/liamcheneyy/Desktop/GCA_001259135.1_2956_6_1_genomic.fna'
# #
# # okay_list = ['A','T','G','C']
# # for filename in glob.iglob('/Users/liamcheneyy/Desktop/untitled folder/*fasta'):
# #     accession = filename.split('/')[-1]
# #     print(accession)
# #     record_dict = {}
# #     for record in SeqIO.parse(filename,'fasta'):
# #         record_dict[record.id] = ""
# #         new_record = ''
# #         for base in record.seq:
# #             if base not in okay_list:
# #                 new_record = new_record + 'N'
# #             else:
# #                 new_record = new_record + base
# #         record_dict[record.id] = new_record
# #
# #     with open('/Users/liamcheneyy/Desktop/fixed/' + accession,'w') as out:
# #         for key in record_dict.keys():
# #             out.write('>' + str(key) + '\n')
# #             out.write(str(record_dict[key]) + '\n')
#
#
# # import glob
# # from shutil import copyfile
# # count = 0
# # for filename in glob.iglob('/Users/liamcheneyy/Desktop/2019-05-27_556_complete_gff/*gff'):
# #     open_list = open('/Users/liamcheneyy/Desktop/Untitled.tsv').read().splitlines()
# #     for line in open_list:
# #         col = line.split('\t')
# #         if filename.split('/')[-1].strip('.gff') == col[0]:
# #             copyfile(filename,'/Users/liamcheneyy/Desktop/gg/' + col[1] + ".gff")
# from time import sleep as sl
#
# keep = []
#
# for i in open('/Users/liamcheneyy/Desktop/xx.txt','r').read().splitlines():
#
#     if '\t' in i:
#         temp_list = []
#         cols = i.split('\t')
#         temp_list.append(cols[0])
#         for i in cols[1:]:
#             new = 'BP1779_' + i
#             temp_list.append(new)
#         fin = '\t'.join(temp_list)
#         keep.append(fin)
#
#     else:
#         if i != 'BP1779_':
#             keep.append(i)
#         else:
#             keep.append('')
#
# for i in keep:
#     print(i)
#

# from ete3 import Tree
#
# # t = Tree('((((H,K)D,(F,I)G)B,E)A,((L,(N,Q)O)J,(P,S)M)C);', format=1)
# t=Tree('((((aaaaaaaaav:1,aaaaaaaaaw:1)1:1,((aaaaaaaaax:1,aaaaaaaaay:1)1:1,(aaaaaaaaaz:1,(aaaaaaaabb:1,(aaaaaaaabc:1,(aaaaaaaabd:1,aaaaaaaabe:1)1:1)1:1)1:1)1:1)1:1)1:1,(aaaaaaaaaa:1,(aaaaaaaaab:1,aaaaaaaaac:1)1:1)1:1)1:1,(((aaaaaaaaad:1,aaaaaaaaae:1)1:1,((aaaaaaaaaf:1,(aaaaaaaaag:1,(aaaaaaaaah:1,(aaaaaaaaai:1,aaaaaaaaaj:1)1:1)1:1)1:1)1:1,(aaaaaaaaak:1,aaaaaaaaal:1)1:1)1:1)1:1,((aaaaaaaaam:1,(aaaaaaaaan:1,(aaaaaaaaao:1,(aaaaaaaaap:1,(aaaaaaaaaq:1,(aaaaaaaaar:1,aaaaaaaaas:1)1:1)1:1)1:1)1:1)1:1)1:1,(aaaaaaaaat:1,aaaaaaaaau:1)1:1)1:1)1:1);', format=1)
#
#
# #creating a trees
# # t=Tree()
# # t.populate(30)
# print(t)
#moving over the trees
# test_root = t.is_root()
# children = t.get_children()[num] #can default is level traversal, can choose either left or right children

#search for strings inside leaf name
# for node in t.traverse():
#         if "av" in node.name:
#             print(node)

# save children results to variables
# ch1,ch2 = t.get_children()

#iterate over leafs of a tree
# for node in t.get_leaves():
#     print(node.name)

#searching the nodes and leaves
# d = t.search_nodes(dist=0.5)[0]
# node_b = t.search_nodes(name="B")[0]
# common_anc = node_b.get_common_ancestor("H","K")

#search leaf by name and get sisters
# leaf_d = t.get_leaves_by_name(name="aaaaaaaaae")[0]
# sisters = leaf_d.get_sisters()[0]

#iterating through leaves
# for leaf in t.iter_leaves(): #returns objects
#     print(leaf.dist)
#
# for leaf in t.iter_leaf_names(): #only returns strings
#     print(leaf)

#adding features to nodes


#writing trees
# print(t.write()) # prints out as a string


# infile = open('/Users/liamcheneyy/Desktop/input.txt','r').read().splitlines()
# add_info = 'SRR8867848_'
#
# save_list = []
# for line in infile:
#     if 'cds' in line:
#         if '\t' not in line:
#             locus = line.split('\t')[0]
#             edit = add_info + locus
#             save_list.append(edit)
#         if '\t' in line:
#             cols = line.split('\t')
#             temp_list = []
#             for j in cols:
#                 j = j.strip('"')
#                 next_edit = (add_info + j)
#                 temp_list.append(next_edit)
#             final = '\t'.join(temp_list)
#             save_list.append(final)
#     else:
#         save_list.append('')
#
# for element in save_list:
#     print(element)

# import pandas as pd
# all_of_interest_path = '/Users/liamcheneyy/Desktop/vcseventh_15/grapetree/all_MGT9_allele_profiles.tsv'
#
# df = pd.read_csv(all_of_interest_path, sep='\t', index_col=False)
#
# allele_count = {}
# for column in df:
#     if '#' not in column:
#         unique_values = list(df[column].unique())
#         alleles_num = len(unique_values)
#         # print(column)
#         # print(unique_values)
#         # sl(1)
#         allele_count[column] = alleles_num
#
# for key, value in allele_count.items():
#     print(str(key) + '\t' + str(value))

# all = open('/Users/liamcheneyy/Desktop/all.txt','r').read().splitlines()
# inp = list(open('/Users/liamcheneyy/Desktop/genes.txt','r').read().splitlines())
#
# for line in all:
#     print(line in inp)
import pandas as pd
from time import sleep as sl

# all = pd.read_csv('/Users/liamcheneyy/Desktop/locus.txt', sep='\t', index_col=False)
# count_dict = all['Locus'].value_counts().to_dict()
#
# for key, value in count_dict.items():
#     print(key, value)

# save_dict = {}
# for key, value in count_dict.items():
#     sub_df = all[all['ST'] == key]
#     strain_num = sub_df.shape[0]
#     if strain_num > 50:
#         sample = sub_df.sample(20)
#         genomes = list(sample['Genome'])
#         save_dict['ST' + str(key)] = genomes
#     else:
#         save_dict['ST' + str(key)] = list(sub_df['Genome'])
#
#
# for key, value in save_dict.items():
#     # print(key)
#     for i in value:
#         print(i)
#     # print()
#     # count.to_csv('/Users/liamcheneyy/Desktop/Untitledoo.csv')
# print(count_dict)

# all = open('/Users/liamcheneyy/Desktop/all.txt', 'r').read().splitlines()
# test = open('/Users/liamcheneyy/Desktop/in.txt', 'r').read().splitlines()
#
# for el in all:
#     if el in test:
#         print("TRUE")
#     else:
#         print("FALSE")

from Bio import SeqIO
# all = open('/Users/liamcheneyy/Desktop/all.txt', 'r').read().splitlines()
#
# save_list = []
# input = SeqIO.parse('/Users/liamcheneyy/Desktop/allele_testinspecies_ref_alleles.fna','fasta')
# for record in input:
#     name = record.id.split(':')[0]
#     if name in all:
#         save_list.append(record)
#
# SeqIO.write(save_list,"/Users/liamcheneyy/Desktop/mgt234_ref_alleles.fna","fasta")

import pandas as pd
# main_df = pd.read_csv('/Users/liamcheneyy/Desktop/complete_metadata.csv', index_col='ID')
# join_df = pd.read_csv('/Users/liamcheneyy/Desktop/simple_blast_results.csv', index_col='ID')
#
# joined = main_df.merge(join_df, left_index=True, right_index=True, how='left')
# joined.to_csv('/Users/liamcheneyy/Desktop/metadata.csv')

# main = open('/Users/liamcheneyy/Downloads/species_MGT.txt', 'r').read().splitlines()
#
# save_list = []
# for line in main:
#     col = line.split('\t')
#     if 'None' not in line:
#         save_list.append(line)
#
# for i in save_list:
#     print(i)

# strains = open('/Users/liamche
from shutil import copyfile

# import glob
# for filename in glob.iglob('/Users/liamcheneyy/Desktop/results/*'):
#     accession = filename.split('/')[-1].split('.')[0]
#     infile = open(filename,'r').read().splitlines()
#     #is tsv
#     tsv = True
#
#     #is csv
#     if len(infile[0].split('\t')) > 1:
#         tsv = True
#
#     #is csv
#     if len(infile[0].split(',')) > 1:
#         tsv = False
#
#     if tsv:
#         with open('/Users/liamcheneyy/Desktop/fixed/' + accession + '_abricate.csv', 'w') as out:
#
#             for element in infile:
#                 col = element.split('\t')
#                 for i in col:
#                     out.write(i + ',')
#                 out.write('\n')
#
#     if not tsv:
#         fixed = accession + '.csv'
#         copyfile(filename, '/Users/liamcheneyy/Desktop/fixed/' + fixed)
#

# infile = open('/Users/liamcheneyy/Desktop/speices_genomes.txt','r').read().splitlines()
# dataset = open('/Users/liamcheneyy/Desktop/MGT_isolate_data.txt').read().splitlines()
#
# inlist = [x for x in infile]
#
# count = 0
# for line in dataset[1:]:
#     col = line.split('\t')
#     acc = col[0]
#     mgt8 = col[7]
#     if acc in inlist:
#         if "None.None" in col[1:7]:

#             print(line)
# for i in infile:
#     acc = i.split(' ')[0]
#     results = i.split('serogroup_cholerae')
#     save_list = []
#     for x in results:
#         if 'rfbV' in x or 'wbfZ' in x:
#             line = x.strip()
#             col = line.split()
#             gene = col[0]
#             percentage = float(col[1])
#             save_list.append([gene,percentage])
#     save_list.sort(key = lambda x:x[1], reverse=True)
#
#     if save_list[0][0] == 'wbfZ_O139':
#         print(acc + '\t' + 'TRUE')
#     if save_list[0][0] == 'rfbV_O1':
#         print(acc + '\t' + 'FALSE')


# infile = open('/Users/liamcheneyy/Desktop/Book5.csv','r').read().splitlines()
#
# for line in infile[1:]:
#     col = line.split(',')
#     strain = col[0]
#
#     ctxB_3 = col[1]
#     ctxB_7 = col[2]
#     ctxB_1 = ctxB_3[3]
#
#     tcpA_el_A226 = col[4]
#     tcpA_el_WT = col[5]
#     tcpA_cla_WT = col[6]
#
#     rstR_cc = col[7]
#     rstR_et = col[8]
#
#     #El Tor
#     # if ctxB_3 == 'TRUE' and tcpA_el_WT == 'TRUE' and rstR_et == 'TRUE':
#     #     print('TRUE')
#
#     #Hatain
    # if ctxB_7 == 'TRUE' and tcpA_el_A226 == 'TRUE' and rstR_et == 'TRUE':
    #     print('TRUE')
#
#     #CLassical
#     # if ctxB_1 == 'FALSE' and tcpA_el_WT == 'FALSE' and rstR_et == 'FALSE':
#     #     print(strain)
#
#     # else:
#     #     print('FALSE')
#
#     # print(col)
#     # sl(1)

# inlist = open('/Users/liamcheneyy/Desktop/arm_abricate.txt').read().splitlines()
#
# inlist = [x for x in inlist]
#
# df = pd.read_csv('/Users/liamcheneyy/Desktop/vcseventh_22/grapetree/seventh/MGT_isolate_data.txt',sep='\t',low_memory=False, index_col=0)
# new_df = df[inlist]
# new_df.to_csv('/Users/liamcheneyy/Desktop/AMR_matrix.csv')
# print(df
#       .shape)
# print(new_df.shape)

# left_df = pd.read_csv('/Users/liamcheneyy/Desktop/left.csv', index_col=0)
# right_df = pd.read_csv('/Users/liamcheneyy/Desktop/Book2.csv', index_col=0)
#
# new = left_df.join(right_df, how='left')
# print(new)
# new.to_csv('/Users/liamcheneyy/Desktop/combined.csv')

# inlist = open('/Users/liamcheneyy/Downloads/genomicepidemiology-resfinder_db-149209df6444/antibiotic_classes.txt','r').read().splitlines()
#
# class_dict = {}
# for i in inlist[1:]:
#     col = i.split('\t')
#     class_dict[col[0]] = []
#     for i in col[1:]:
#         if i != '':
#             class_dict[col[0]].append(i)

# inlist = open('/Users/liamcheneyy/Desktop/Book3.txt','r').read().splitlines()
# inlist_p = [x.lower() for x in inlist]
#
# for line in open('/Users/liamcheneyy/Desktop/amr/card_database.txt').read().splitlines():
#     col = line.split('\t')
#     gene = col[0]
#     Class  = col[3].replace('antibiotic','').lower()
#
#     if Class not in inlist_p:
#         print(Class)
#         print(inlist_p)

# infile = open('/Users/liamcheneyy/Desktop/results.txt').read().splitlines()
#
# for line in infile:
#     strain = line.split(' ')[0]
#     count = line.count('biotype_cholerae')
#     if count > 1:
#         cols = line.split('biotype_cholerae')[1:]
#         save = []
#         for i in cols:
#             # print(cols)
#             if 'GCA' not in strain:
#                 keep = i.split('contig')[0].strip()
#                 save.append(keep)
#             if 'GCA' in strain:
#                 keep = i.split('NZ_')[0].strip()
#                 save.append(keep)
#
#         save_string = '-'.join(save)
#         print(strain, ',', save_string)
#
#     else:
#         if len(line.split(' ')) > 1:
#             single_keep = line.split('biotype_cholerae')[1].split('contig')[0].strip()
#             print(strain, ',', single_keep)
#         if 'biotype' not in line:
#             print(strain)

# infile = open('/Users/liamcheneyy/Desktop/results.csv').read().splitlines()
# infile_dict = {}
# for line in infile:
#     col = line.split(',')
#     infile_dict[col[0]] = col[1:]
#
# df = pd.read_csv('/Users/liamcheneyy/Desktop/combinations.csv', index_col=0)
#
# save_dict = {}
# for column in df:
#     save_dict[column] = {}
#     sub = df[df[column] == True]
#     strains = list(sub.index)
#     for i in strains:
#         save_dict[column][i] = infile_dict[i]
#
# for key, value in save_dict.items():
#     with open('/Users/liamcheneyy/Desktop/combs/' + str(len(value)) + '_' + key + '.txt','w') as out:
#         for i in value:
#             out.write(i + '\t')
#             for x in value[i]:
#                 out.write(x + '\t')
#             out.write('\n')

# import glob
#
# for filename in glob.iglob('/Users/liamcheneyy/Desktop/combs/1_*txt'):
#     infile = open(filename,'r').read().splitlines()
#     combination = '_'.join(filename.split('_')[1:]).replace('.txt','')
#     for line in infile:
#         print(combination, line)

# infile = open('/Users/liamcheneyy/Desktop/combs/all.txt').read().splitlines()
#
# for line in infile:
#     col = line.split('\t')
#     rstR = col[-1]
#     stain = col[0]
#
#     if rstR.count('100') == 2:
#         print(stain,sep='\t')
#
#     if '-' in rstR:
#         if rstR.count('100') != 2:
#             print(stain,sep='\t')

# import glob
#
# save_dict = {}
# for filename in glob.iglob('/Users/liamcheneyy/Desktop/out/*txt'):
#     infile = open(filename,'r').read().splitlines()
#     strain = filename.split('/')[-1].strip('.txt')
#     keep = infile[0].split('virulence_cholerae')[0]
#     keep = keep.split('biotype_cholerae')[1:]
#     save_list = []
#     save_percen = []
#     save_cov = []
#     for i in keep:
#         ctxb = i[:23].strip(' ')
#         save_list.append(ctxb)
#         percen = ctxb.split()[1]
#         save_percen.append(percen)
#         cov = ctxb.split()[2]
#         save_cov.append(cov)
#
#     if save_percen[0] == save_percen[1] and (save_cov[0] == '375'):
#         # print(strain, save_percen, save_cov)
#         print(strain)
#
#     save_dict[strain] = '-'.join(save_list)
#     # for i in

# for line in open('/Users/liamcheneyy/Desktop/Book4.csv').read().splitlines():
#
#     col = line.split(',')
#     strain = col[0]
#     classes = col[3].count(';')
#     print(classes)

    # gene = col[0]
    # Class  = col[3].replace('antibiotic','').lower()
    #
    # if Class not in inlist_p:
    #     print(Class)
    #     print(inlist_p)

 # import pandas as pd
# df = pd.read_csv('/Users/liamcheneyy/Desktop/combinations.csv')
#
# save_dict = {}
# for column in df:
#     sub = df[df[column] == True]
#     strains = list(sub.index)
#
#     print(sub.shape)

# import glob,re
#
# save_dict = {}
# for file in glob.iglob('/Users/liamcheneyy/Desktop/x/*txt'):
#     infile = open(file).read()
#     sep = re.split("\n[0-9]*\:", infile)
#
#     bioproject = file.split('/')[-1].strip('.txt')
#
#     for strains in sep:
#         strain = ''
#         year = ''
#         country = ''
#         strain_id = ''
#
#         lines = strains.split('\n')
#
#         for line in lines:
#             if 'SRA:' in line:
#                 strain = line.split('SRA: ')[-1]
#
#                 if strain not in save_dict.keys():
#                     save_dict[strain] = {}
#
#                 elif strain in save_dict.keys():
#                     print(strain)
#
#             if ('date' in line or 'year' in line) and ('update' not in line) and ('issing' not in line) and ('not' not in line):
#                 year = line.split('=')[-1].strip(""" " """)
#
#             if 'location' in line and ('long' not in line and 'lat' not in line and 'issing' not in line) and ('not' not in line):
#                 if 'not' in line:
#                     pass
#                 else:
#                     country = line.split('=')[-1].strip(""" " """)
#
#             if 'strain' in line and ('not' not in line) and ('issing' not in line):
#                 strain_id = line.split('=')[-1].strip(""" " """)
#
#             # if 'isolate' in line:
#             #     print(line)
#
#         save_dict[strain] = {'country':country, 'year':year, 'strain':strain_id, 'bioproject':bioproject}
#
# import pandas as pd
#
# df = pd.DataFrame(save_dict).T
# df.replace()
#
# df.to_csv('/Users/liamcheneyy/Desktop/x/all.csv')

# keep_list = []
# df = pd.read_csv('/Users/liamcheneyy/Desktop/Untitle.txt', sep='\t', index_col=0)
#
# count = 1
#
# for col in df:
#     sub = df[df[col] == True]
#     sub = sub[col]
#     strains = list(sub.index)
#     if '/' in col:
#         col = col.split('/')[0] + '-' + col.split('/')[-1]
#
#     if len(strains) <= 21:
#         keep_list.append(strains)
#
#     # sub.to_csv('/Users/liamcheneyy/Desktop/combs/' + str(len(strains)) + '_' + str(col) + '.csv', header=['ID',col])
#
# keep = [i for x in keep_list for i in x]
# for i in keep:
#     print(i)

# import glob
#
# details = pd.read_csv('/Users/liamcheneyy/Desktop/test/simple_blast_results.csv', index_col=0)
#
# for filename in glob.iglob('/Users/liamcheneyy/Desktop/combs/*csv'):
#     file = filename.strip('.csv')
#     df = pd.read_csv(filename, index_col=0)
#     strains = list(df.index)
#     sub = details[details.index.isin(strains)]
#     sub = sub.replace(to_replace='NaN',value='')
#     sub.to_csv(file + '_details.csv')
#     print(sub)

# import glob
# for filename in glob.iglob('/Users/liamcheneyy/Desktop/taking/*'):
#     infile = open(filename).read().splitlines()
#     strain = filename.split('/')[-1].strip('.tsv')
#     save_list = []
#     for line in infile:
#         if 'biotype' in line and 'rstR' in line:
#             col = line.split('\t')
#             cov = col[3]
#             percen = int(cov.split('/')[0].strip('')) / int(cov.split('/')[0].strip(''))
#             save_list.append(percen)
#
#     if save_list[0] == save_list[1]:
#         print(strain)

from Bio import SeqIO
from Bio.Seq import Seq


for record in SeqIO.parse('/Users/liamcheneyy/Desktop/tt/Untitled.fasta','fasta'):
    outfile = '/Users/liamcheneyy/Desktop/ref_index/xx.fasta'
    SeqIO.write(record,outfile,'fasta')


