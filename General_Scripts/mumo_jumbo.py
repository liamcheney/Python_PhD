import glob
from time import sleep as sl
import matplotlib
# matplotlib.use("TkAgg")
# from matplotlib import pyplot as plt
# from collections import Counter
#
# infile_path = '/Users/liamcheneyy/Desktop/failure_reasons/'

# genome_size_list = []
# large_contig_list = []
# kraken_list = []
# n50_list = []
#
# contam_count = 0
# gen_count = 0
# for filename in glob.iglob(infile_path + '*'):
#     infile = open(filename).read().splitlines()
#     file_acc = filename.split('/')[-1].strip('_failure_reason.txt')
#     for line in infile:
#         if "Genome outside of allowed range:" in line:
#             genome_size = line.split('(')[-1].strip(')')
#             genome_size_list.append([file_acc,float(genome_size)])
#             # gen_count = gen_count + 1
#
#         if "largest contig <" in line:
#             l_contig = line.split('(')[-1].strip(')')
#             large_contig_list.append(int(l_contig))
#             gen_count = gen_count + 1
#
#
#         if "N50 less" in line:
#             n50 = int(line.split('(')[-1].strip(')')) * 100
#             n50_list.append(int(n50))
#             # gen_count = gen_count + 1
#
#
#         if "contaminated" in line:
#             contam_species = line.split('\t')[-1]
#             kraken_list.append(contam_species)
#             contam_count = contam_count + 1
#
# for i,b in genome_size_list:
#
#     print(i)
# plt.hist(kraken_list, bins=40)
# plt.xticks(rotation=90)
# plt.show()
#
# kraken = Counter(kraken_list)
#
# import plotly.plotly as py
# import plotly.graph_objs as go
#
# import pandas as pd
#
# df = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/2014_us_cities.csv')
#
# df['text'] = df['name'] + '<br>Population ' + (df['pop']/1e6).astype(str)+' million'
# limits = [(0,2),(3,10),(11,20),(21,50),(50,3000)]
# colors = ["rgb(0,116,217)","rgb(255,65,54)","rgb(133,20,75)","rgb(255,133,27)","lightgrey"]
# cities = []
# scale = 5000
# #
# # for i in range(len(limits)):
# #     lim = limits[i]
# #     df_sub = df[lim[0]:lim[1]]
# #     city = go.Scattergeo(
# #         locationmode = 'USA-states',
# #         lon = df_sub['lon'],
# #         lat = df_sub['lat'],
# #         text = df_sub['text'],
# #         marker = go.scattergeo.Marker(
# #             size = df_sub['pop']/scale,
# #             color = colors[i],
# #             line = go.scattergeo.marker.Line(
# #                 width=0.5, color='rgb(40,40,40)'
# #             ),
# #             sizemode = 'area'
# #         ),
# #         name = '{0} - {1}'.format(lim[0],lim[1]) )
# #     cities.append(city)
# #
# # layout = go.Layout(
# #         title = go.layout.Title(
# #             text = '2014 US city populations<br>(Click legend to toggle traces)'
# #         ),
# #         showlegend = True,
# #         geo = go.layout.Geo(
# #             scope = 'usa',
# #             projection = go.layout.geo.Projection(
# #                 type='albers usa'
# #             ),
# #             showland = True,
# #             landcolor = 'rgb(217, 217, 217)',
# #             subunitwidth=1,
# #             countrywidth=1,
# #             subunitcolor="rgb(255, 255, 255)",
# #             countrycolor="rgb(255, 255, 255)"
# #         )
# #     )
# #
# # fig = go.Figure(data=cities, layout=layout)
# # py.iplot(fig, filename='d3-bubble-map-populations')

# infile_path = open('/Users/liamcheneyy/Desktop/Untitled.txt').read().splitlines()
#
# issue_list = ["virulence","chemotaxis","biotic","pilus",'membrane','secretion','t6ss','t1ss','t2ss','toxin','periplasmic','toxin','periplasmic','extracellular','capsid','capsular','wall'
#               ,'flagella','mobile','plasmid','pilin','drug','receptor','resistance','peptidoglycan','transposase','twin','fimbrial','transporter']
#
# non_issue_dict = {}
# issue_dict = {}
# for line in infile_path:
#     cell = line.split('\t')
#     for el in issue_list:
#         if el in cell[1].lower():
#             issue_dict[cell[0]] = cell[1]
#         else:
#             non_issue_dict[cell[0]] = cell[1]
#
# # print(len(issue_dict.keys()))
# for key, value in issue_dict.items():
#     print(key + '\t' + value)

# infile = open('/Users/liamcheneyy/Desktop/vc_core_genes.txt','r').read().splitlines()
# bio_list= list([x for x in open('/Users/liamcheneyy/Desktop/biocycle_list.txt').read().splitlines()])
#
# for line in infile:
#     if line in bio_list:
#         print(line + '\t' + "TRUE")
#     else:
#         print(line + '\t' + "FALSE")

bio_list= list([x for x in open('/Users/liamcheneyy/Desktop/biocycle_list.txt').read().splitlines()])
