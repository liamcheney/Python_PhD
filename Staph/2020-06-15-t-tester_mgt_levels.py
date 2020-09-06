##reading in files
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from pylab import figure, text, scatter, show
import numpy as np
import matplotlib.mlab as mlab
from time import sleep as sl

def get_mgt9_values(info, attr):
    mgt9_values = list(info[attr])
    mgt9_values = [x for x in mgt9_values if (isinstance(x, float) or isinstance(x, int))]
    return mgt9_values

def mean_confidence_interval(data, confidence=0.95):

    size = len(data)
    t = stats.t(df = size - 1).ppf((1 + confidence) /2)
    error = t * np.std(data, ddof=1) / np.sqrt(size)
    means = np.mean(data)
    return means, error

def t_tester(mgt9_values, values, mgt_level):
    # tStat, pValue = stats.ttest_ind(values, mgt9_values, equal_var=False)  # run independent sample T-Test
    # print(f"{mgt_level} P-Value:{pValue}")

    tStat, pValue = stats.mannwhitneyu(values,mgt9_values)
    # print(f"{mgt_level} P-Value:{pValue}")

    pValue = round(pValue,6)
    return pValue

def mgt_level_iter(loci_dict, info, mgt9_values):
    pValue_dict = {}
    CI_dict = {}
    all_values_list = []
    for mgt_level in loci_dict.keys():
        if mgt_level != 'mgt9':
            loci = loci_dict[mgt_level]
            sub_df = info[info[attr].index.isin(loci)]
            values_pre = list(sub_df[attr])
            values = []
            for x in values_pre:
                if not isinstance(x, str):
                    values.append(x)

            all_values_list.append(values)
            # print(mgt_level, len(values_pre), len(values))

            ##Confidence intervals
            mean, error = mean_confidence_interval(values,0.95)
            CI_dict[mgt_level] = [mean, error]
            # print(mgt_level, *values)
            # print(mgt_level, mean, error, mean-error, mean+error)

            ##t-testing
            # pValue = t_tester(mgt9_values, values, mgt_level)
            # pValue_dict[mgt_level] = pValue

    return pValue_dict, CI_dict, all_values_list

def plotting(pValue_dict, CI_dict, all_values_list,attr):

    #########P-values
    fig1, ax1 = plt.subplots()
    ax1.set_title(f"Comparing MGT level differences for {attr} ")
    ax1.boxplot(all_values_list, showmeans=True)

    locs, labels = plt.xticks()
    plt.xticks(locs, ('MGT3', 'MGT4', 'MGT5', 'MGT6', 'MGT7','MGT8','MGT9'))
    # plt.ylim(0,40)

    if pValue_dict:
        x_point = 0.7
        for el in pValue_dict.values():
            text(x_point, .5, f"{el}", fontsize=8)
            x_point += 1

    plt.savefig('/Users/liamcheneyy/Desktop/dists/' + attr + '_boxplot_pvalues.png')
    plt.clf()

    # #########Convidence Intervals
    means_list = [x[0] for x in CI_dict.values()]
    error_list = [x[1] for x in CI_dict.values()]
    levels = [x for x in CI_dict.keys()]
    print(levels)
    print(means_list)
    print(error_list)
    fig, ax = plt.subplots()
    ax.set_title(f'{attr} Confidence Interval')
    plt.errorbar(x=means_list, y=levels, xerr=error_list, fmt='o', color='k')
    mgt9_mean = np.mean(all_values_list[-1])
    # print(mgt9_mean)
    # print(len(all_values_list[-1]))
    plt.axvline(mgt9_mean, ls='--')
    plt.savefig('/Users/liamcheneyy/Desktop/dists/' + attr + '_CI.png')
    plt.clf()

    ##########Hists
    titles = ['MGT3', 'MGT4', 'MGT5', 'MGT6', 'MGT7','MGT8','MGT9']
    f, a = plt.subplots(2, 3)
    a = a.ravel()

    for idx, ax in enumerate(a):
        kwargs = dict(bins=30)
        ax.hist(all_values_list[idx], **kwargs, range=(sorted(all_values_list[-2])[0],sorted(all_values_list[-2])[-1]))
        ax.set_title(titles[idx])

    plt.tight_layout()
    plt.savefig('/Users/liamcheneyy/Desktop/dists/' + attr + '_hists.png')

def outstats(all_values_list):
    cvar = 3
    for i in all_values_list:
        MGT = 'MGT' + str(cvar)
        print(MGT, np.mean(i))
        with open('/Users/liamcheneyy/Desktop/loci/attr_' + str(MGT) + '.txt', 'w') as out:
            out.write(str(MGT) + '\n')
            for x in i:
                out.write(str(x) + '\n')
        cvar += 1



if __name__ == '__main__':
    ###inputs
    mgt3 = open('/Users/liamcheneyy/Desktop/loci/MGT3_gene_accessions.txt').read().splitlines()
    mgt4 = open('/Users/liamcheneyy/Desktop/loci/MGT4_gene_accessions.txt').read().splitlines()
    mgt5 = open('/Users/liamcheneyy/Desktop/loci/MGT5_gene_accessions.txt').read().splitlines()
    mgt6 = open('/Users/liamcheneyy/Desktop/loci/MGT6_gene_accessions.txt').read().splitlines()
    mgt7 = open('/Users/liamcheneyy/Desktop/loci/MGT7_gene_accessions.txt').read().splitlines()
    mgt8 = open('/Users/liamcheneyy/Desktop/loci/MGT8_gene_accessions.txt').read().splitlines()
    mgt9 = open('/Users/liamcheneyy/Desktop/loci/MGT9_gene_accessions.txt').read().splitlines()

    input_path = "/Users/liamcheneyy/Desktop/Preferences.xlsx"
    info = pd.read_excel(open(input_path, 'rb'), index_col=0, sheet_name='Preferences').fillna("none")

    ##variables
    loci_dict = {'mgt3':mgt3,'mgt4':mgt4,'mgt5':mgt5,'mgt6':mgt6,'mgt7':mgt7,'mgt8':mgt8,'mgt9':mgt9}
    # to_do_list = ['log_dNdS', 'Recombination Events']
    to_do_list = ['Allele_Changes_Rate']

    for i in to_do_list:
        attr = i

        mgt9_values = get_mgt9_values(info, attr)

        pValue_dict, CI_dict, all_values_list = mgt_level_iter(loci_dict, info, mgt9_values)
        all_values_list.append(mgt9_values)
        plotting(pValue_dict, CI_dict, all_values_list, attr)

        outstats(all_values_list)