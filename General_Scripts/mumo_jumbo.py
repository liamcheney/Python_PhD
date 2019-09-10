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


infile = open('/Users/liamcheneyy/Desktop/input.txt','r').read().splitlines()
add_info = 'SRR8867848_'

save_list = []
for line in infile:
    if 'cds' in line:
        if '\t' not in line:
            locus = line.split('\t')[0]
            edit = add_info + locus
            save_list.append(edit)
        if '\t' in line:
            cols = line.split('\t')
            temp_list = []
            for j in cols:
                j = j.strip('"')
                next_edit = (add_info + j)
                temp_list.append(next_edit)
            final = '\t'.join(temp_list)
            save_list.append(final)
    else:
        save_list.append('')

for element in save_list:
    print(element)

import pandas as pd
all_of_interest_path = '/Users/liamcheneyy/Desktop/vcseventh_15/grapetree/all_MGT9_allele_profiles.tsv'

df = pd.read_csv(all_of_interest_path, sep='\t', index_col=False)

allele_count = {}
for column in df:
    if '#' not in column:
        unique_values = list(df[column].unique())
        alleles_num = len(unique_values)
        # print(column)
        # print(unique_values)
        # sl(1)
        allele_count[column] = alleles_num

for key, value in allele_count.items():
    print(str(key) + '\t' + str(value))