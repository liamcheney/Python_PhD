import pandas as pd
from time import sleep as sl
from collections import Counter

#stripts out cds, vc and gene info from gff file
# wrote to file as useful in other analysis, also takes long time to continually create
gff = open('/Users/liam/Desktop/reference/GCA_000006745.gff', 'r').read().splitlines()
gene_cds_locus = open('/Users/liam/Desktop/gene_cds_locus.csv','w')
for i in gff:
    if 'ID=gene' in i and 'RefSeq' in i:
        gencol = i.split('\t')
        gene = gencol[8].split(';')[0].split('_')[0][3:]
        locus = gencol[8].split(';')[-1][10:]
        for j in gff:
            if 'ID=cds' in j:
                cdscol = j.split('\t')
                start = cdscol[3]
                end = cdscol[4]
                strand = cdscol[6]
                pgene = cdscol[8].split(';')[1].split('=')[1]
                cds = cdscol[8].split(';')[0].split('_')[0][3:]
                if gene == pgene:
                    gene_cds_locus.write(gene + ',' + cds + ',' + locus + ',' + start + ',' + end + ',' + strand + '\n')
                    print(cdscol)
#
# ###reads in an excel file and selects a list of core genes for D3
# cg_list = pd.read_excel('/Users/liam/Desktop/All_Core_Genome/alldatasets_coregenome.xlsx', None)
# d3_cg = cg_list['D3 conversion'].iloc[:,6].dropna()
# d3_cg_list = d3_cg.values.tolist()
# total_cg = len(d3_cg_list)
#
# ### as patric does not have the VC, annotation and subsystem in same table, script extracts extracts the system and inserts
# all_genome_features = pd.read_csv('/Users/liam/Desktop/patric/PATRIC_genome_feature.csv')
# all_subsystems = pd.read_csv('/Users/liam/Desktop/patric/PATRIC_all_subsystems.csv')
# merged = all_subsystems.merge(all_genome_features, left_on = 'PATRIC ID', right_on = 'PATRIC ID')
#
# d3_cds_sub = merged[merged['RefSeq Locus Tag'].isin(d3_cg_list)]
# # d3_cds_sub.to_csv('/Users/liam/Desktop/patric/PATRIC_d3_cds_sub.csv', index=False)
#
# #create list of subclass types and calculate percentage each sub class is for core genes
# d3_subclass_list = d3_cds_sub.iloc[:,0].tolist()
# subsys_type = set(d3_subclass_list)
# print(len(d3_subclass_list))
# d3_subsytems = open('/Users/liam/Desktop/patric/PATRIC_d3_subsystems.csv','w')
# un_annotated = 0
# for i in subsys_type:
#     percen = int((d3_subclass_list.count(i)))
#     un_annotated = un_annotated + percen
#     d3_subsytems.write(str(i) + ',' + str(percen) + '\n')
# d3_subsytems.write('Un-annotated' + ',' + str(un_annotated) + '\n')