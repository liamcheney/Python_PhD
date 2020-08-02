##reading in files
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import matplotlib.pyplot as plt
from pylab import figure, text, scatter, show


mgt2 = open('/Users/liamcheney/Desktop/loci/MGT2_gene_accessions.txt').read().splitlines()
mgt3 = open('/Users/liamcheney/Desktop/loci/MGT3_gene_accessions.txt').read().splitlines()
mgt4 = open('/Users/liamcheney/Desktop/loci/MGT4_gene_accessions.txt').read().splitlines()
mgt5 = open('/Users/liamcheney/Desktop/loci/MGT5_gene_accessions.txt').read().splitlines()
mgt6 = open('/Users/liamcheney/Desktop/loci/MGT6_gene_accessions.txt').read().splitlines()
mgt7 = open('/Users/liamcheney/Desktop/loci/MGT7_gene_accessions.txt').read().splitlines()
mgt8 = open('/Users/liamcheney/Desktop/loci/MGT8_gene_accessions.txt').read().splitlines()

input_path = "/Users/liamcheney/Desktop/Preferences.xlsx"
info = pd.read_excel(open(input_path, 'rb'), index_col=0, sheet_name='Preferences').fillna("none")

loci_dict = {'mgt2':mgt2,'mgt3':mgt3,'mgt4':mgt4,'mgt5':mgt5,'mgt6':mgt6,'mgt7':mgt7,'mgt8':mgt8}

# attr = 'dNdS'
# attr = 'Length'
# attr = 'Allelic_Changes'
attr = 'Recombination Events'

attri_df = info[attr]
mgt8_values = list(info[attr])
mgt8_values = [x for x in mgt8_values if isinstance(x,float) or isinstance(x,int)]

save_list = []
pvalue_list = []
for mgt_level in loci_dict.keys():
    lists = []
    if mgt_level != 'mgt8':
        loci = loci_dict[mgt_level]
        sub_df = info[info[attr].index.isin(loci)]
        values = list(sub_df[attr])
        for i in values:
            if i != 'none':
                lists.append(i)

        # if mgt_level == 'mgt5':
        #     for i in values:
        #         print(i)

        ##runing t-test
        # tStat, pValue = stats.ttest_ind(lists, mgt8_values, equal_var=False)  # run independent sample T-Test
        # print(f"{# mgt_level} P-Value:{pValue}")

        tStat, pValue = stats.mannwhitneyu(lists,mgt8_values)
        print(f"{mgt_level} P-Value:{pValue}")

        pValue = round(pValue,6)
        save_list.append(lists)
        pvalue_list.append(pValue)

        with open('/Users/liamcheney/Desktop/loci/' + attr + '_' + mgt_level + '.txt','w') as out:
            for i in lists:
                out.write(str(i) + '\n')
#
save_list.append(mgt8_values)
data = save_list
fig1, ax1 = plt.subplots()
ax1.set_title(f"Comparing MGT level differences for {attr} ")
ax1.boxplot(data)

locs, labels = plt.xticks()
plt.xticks(locs, ('MGT2', 'MGT3', 'MGT4', 'MGT5', 'MGT6','MGT7','MGT8'))
plt.ylim(0,40)

stat = 0.7
for el in pvalue_list:
    text(stat, 37 , f"{el}", fontsize=8)
    stat += 1

plt.show()

