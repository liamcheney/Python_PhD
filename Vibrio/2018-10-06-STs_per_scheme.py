import pandas as pd
from time import sleep as sl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pandas import ExcelWriter

import seaborn as sns
sns.set_palette(sns.color_palette("tab20", 20))

from make_outp_excel_ready import split_hgts
mpl.rcParams['pdf.fonttype'] = 42


def data_metadata_import_merge(datafile,metafile,precombined):
    if precombined != "":
        fulldata = pd.read_excel(precombined, index_col=0)
    else:
        mgt_data = split_hgts(datafile, False)
        metadata = pd.read_excel(metafile,index_col=0, dtype=str,na_values="0")
        fulldata = pd.merge(mgt_data,metadata,left_index=True,right_index=True)

    return fulldata

# def major_st_over_time(dataframe):
#     df = pd.read_csv("/Users/michaelpayne/Documents/UNSW/OneDrive - UNSW/HGT_Paper/figures/dt104_major_ST_over_time/major_st_perc_over_time_cc.txt",sep="\t",header=0).set_index('Year')
#     df = df[4:]
#
#     plot = df.plot.line(figsize=(10,6))
#
#     # ax.legend(loc=1)
#     # plt.ylabel("Strain count")
#     # plt.xticks(df.index)
#     plot.set_yticks([0,20,40,60,80,100])
#     plot.set_xticks(range(1990,2011,2))
#     plot.set_xlim(1990,2018)
#     plot.legend(loc=0)
#     #plt.show
#
#     plot = plot.get_figure()
#     # plot.tight_layout()
#     plot.show()

def change_recursive(change_dict,string,count,no):
    if count == no:
        return string
    else:
        l = string.split("-")
        tocheck = l[count]
        tochange = count-1
        l = l[:tochange] + [change_dict[count][tocheck]] + l[tochange+1:]
        nstring = "-".join(l)
        count += 1
        return change_recursive(change_dict,nstring,count,no)


def get_inconsistent(all,scheme):

    nos = all.groupby(scheme).count()
    nos = (nos["Strain"])
    no = int(scheme[-1])

    ind_lis = list(nos.T.index)

    # print(ind_lis)
    exist = {}
    exist2 = {}
    for pos in range(1,no):
        exist[pos] = {}
        exist2[pos] = {}
        for i in ind_lis:
            l = i.split("-")
            id = l[pos]
            up = l[pos-1]
            if id not in exist[pos]:
                exist[pos][id] = [up]
            else:
                exist[pos][id].append(up)
        for i in exist[pos]:
            idlis = exist[pos][i]
            nums1 = list(set(exist[pos][i]))
            nums = [(x,idlis.count(x)) for x in nums1]
            mx = max(nums, key=lambda item: item[1])

            # exist2[pos][i] = mx[0]
            if len(nums1) > 1:
                # exist2[pos][i] = mx[0]+"*"
                exist2[pos][i] = "/".join(nums1)
            # elif len(nums1) > 1 and len(nums1) <= 3 :
            #     exist2[pos][i] = "/".join(nums1)
            else:
                exist2[pos][i] = mx[0]

    # print(exist2[1])
    conv = dict()
    convr = {}
    ind_lis2 = list(nos.T.index)
    for i in ind_lis2:
        newi = change_recursive(exist2,i,1,no)
        convr[i] = newi
        if newi in conv:
            conv[newi].append(i)
        else:
            conv[newi] = [i]
    # print(convr)
    return convr





def leekit_country(leekit,schemes,outfolder):
    for scheme in schemes:

        scheme = "id_string_"+scheme

        print(leekit.columns)
        df2 = leekit.groupby(["Location (Corrected - Country)", scheme])["Location (Corrected - Country)"].count().unstack(scheme).fillna(0)
        print(df2.columns)
        df2 = df2[df2.sum().sort_values(ascending=False).index]
        print(df2.columns)
        df2 = df2.reindex(["Austria","Germany","Denmark","Netherlands","Luxembourg","Switzerland","France","Spain","Ireland","Poland","Czech Republic","Israel","Morocco","Thailand","Taiwan","New Zealand","Canada","United States","Argentina"])
        print(df2.columns)
        # get_inconsistent(df2)

        # print(list(df2.T.index))

        ###### with singletons separate plot

        # plt = df2.plot.bar(stacked=True,figsize=(10,8),width=0.8).get_figure()
        # # plt.tight_layout()
        # plt.savefig(scheme+"_leekit_country_distribution.pdf")
        # plt.show()

        #####################

        ###### with singletons together plot
        # df2.to_excel(writer, 'leekit_country_' + scheme)
        series = df2.sum(axis=0)
        mask = (series).gt(1)

        # print(series)
        tokeep = list(series[mask].index)
        tocombine = list(series[~mask].index)

        morethanone = df2[tokeep]

        singles = pd.DataFrame(columns=["Singletons"])

        singles["Singletons"] = df2[tocombine].sum(axis=1)

        final = pd.concat([morethanone,singles],axis=1)







        plt2 = final.plot.bar(stacked=True,figsize=(10,8),width=0.8).get_figure()
        # plt2.xlabel("Country")
        # plt2.ylabel("Strain count")
        # plt2.tight_layout()

        plt2.savefig(outfolder+scheme+"_leekit_country_distribution_no1.pdf")

        ##################
def all_source(total,schemes,outfolder):
    schemes+=["MGT92","MGT95","MGT910"]
    for scheme in schemes:
        if len(scheme) == 4:
            scheme = "id_string_"+scheme

        to_keep = ["clinical", "environmental/other"]
        total = total[total["met2"].isin(to_keep)]

        to_keep = ["United Kingdom: Scotland", "United Kingdom: England", "United Kingdom: Wales"]
        total = total[total["Location (Corrected - Country)"].isin(to_keep)]

        df2 = total.groupby(['met2', scheme])['met2'].count().unstack(scheme).fillna(0)

        df2 = df2[df2.sum().sort_values(ascending=False).index]

        # df2 = df2.reindex(["Austria","Germany","Denmark","Netherlands","Luxembourg","Switzerland","France","Spain","Ireland","Poland","Czech Republic","Israel","Morocco","Thailand","Taiwan","New Zealand","Canada","United States","Argentina"])

        # get_inconsistent(df2)

        # print(list(df2.T.index))

        ###### with singletons separate plot

        # plt = df2.plot.bar(stacked=True,figsize=(10,8),width=0.8).get_figure()
        # # plt.tight_layout()
        # plt.savefig(scheme+"_leekit_country_distribution.pdf")
        # plt.show()

        #####################

        ###### with singletons together plot
        # df2.to_excel(writer, 'leekit_country_' + scheme)
        series = df2.sum(axis=0)
        mask = (series).gt(1)

        # print(series)
        tokeep = list(series[mask].index)
        tocombine = list(series[~mask].index)

        morethanone = df2[tokeep]

        singles = pd.DataFrame(columns=["Singletons"])

        singles["Singletons"] = df2[tocombine].sum(axis=1)

        final = pd.concat([morethanone,singles],axis=1)







        plt2 = final.plot.bar(stacked=True,figsize=(10,8),width=0.8).get_figure()
        # plt2.xlabel("Country")
        # plt2.ylabel("Strain count")
        # plt2.tight_layout()

        plt2.savefig(outfolder+scheme+"_with_string_source_distribution_no1.pdf")

        ##################

def mather_UK_time(mather,schemes,outfolder):
    # print(mather)
    to_keep = ["United Kingdom: Scotland","United Kingdom: England","United Kingdom: Wales"]
    mather = mather[mather["Location (Corrected - Country)"].isin(to_keep)]
    # print(mather)
    for scheme in schemes:
        scheme = "id_string_" + scheme
        df2 = mather.groupby(['Collection Year', scheme])['Collection Year'].count().unstack(scheme).fillna(0)
        # print(df2)

        # df2 = df2.reindex(df2.T.sum().sort_values().T.index,axis=1)
        df2 = df2[df2.sum().sort_values(ascending=False).index]
        # df2.to_excel(writer, 'mather_time_' + scheme)
        series = df2.sum(axis=0)
        mask = (series).gt(1)

        # print(series)
        tokeep = list(series[mask].index)
        tocombine = list(series[~mask].index)

        morethanone = df2[tokeep]

        singles = pd.DataFrame(columns=["Singletons"])

        singles["Singletons"] = df2[tocombine].sum(axis=1)

        final = pd.concat([morethanone, singles], axis=1)






        plt2 = final.plot.area(stacked=True, figsize=(8, 6),lw=0,fontsize=12).get_figure()
        # plt2.xlabel("Year")
        # plt2.ylabel("Strain count")
        plt2.tight_layout()

        plt2.savefig(outfolder+scheme + "_mather_time_UK_distribution_no1.pdf")

def major_ST_over_time(all,schemes,outfolder,name,cc_or_st):
    df = pd.DataFrame()
    c = 0
    for scheme in schemes:
        scheme = "id_string_" + scheme
        # all[['Collection Year']].apply(pd.to_numeric, errors='ignore')

        df2 = all.groupby(['Collection Year', scheme])['Collection Year'].count().unstack(scheme).fillna(0)
        # print(df2)
        # df3 = list(df2.sum(axis=0).sort_values(ascending=False).T.index)
        maxst = list(df2.sum().sort_values(ascending=False).index)[0]
        df2 = df2.div(df2.sum(axis=1), axis=0).multiply(100)



        # print(df2)
        # sl(10)

        if c == 0:
            df = df2[maxst]
        else:
            df = pd.concat([df, df2[maxst]], axis=1)


        c+=1
        # df2 = df2.reindex(df3)
    # df = df.drop(["??"])

    # print(df)

    # df.to_csv("/Users/michaelpayne/Documents/UNSW/OneDrive - UNSW/HGT_Paper/figures/dt104_major_ST_over_time/3-major_cc_perc_over_time_rerun2.txt",sep="\t")


    # df.to_excel(writer, 'major_over_time')


    # df = df.drop(df.columns[[7,8]], axis=1)
    # print(df)

    # print(df)

    plot = df.plot.line(figsize=(10, 6))

    # ax.legend(loc=1)
    # plt.ylabel("Strain count")
    # plt.xticks(df.index)
    plot.set_yticks([0, 20, 40, 60, 80, 100])
    if name == "DT104":
        plot.set_xticks(range(1990, 2011, 2))
        plot.set_xlim(1990, 2018)
    elif name == "DT160":
        plot.set_xticks(range(1998, 2012, 2))
        plot.set_xlim(1998, 2012)
    plot.legend(loc=0)
    plot = plot.get_figure()
    plot.tight_layout()
    plot.savefig(outfolder+ cc_or_st + "_" +
        name +"_major_st_perc_over_time_rerun.pdf")

        # print(df2)

def add_singletones(df): ####creating for graph, just a large block for the single STs instead of lots of little blocks.
    df = df[df.sum().sort_values(ascending=False).index]
    series = df.sum(axis=0)
    mask = (series).gt(1)

    # print(series)
    tokeep = list(series[mask].index)
    tocombine = list(series[~mask].index)

    morethanone = df[tokeep]

    singles = pd.DataFrame(columns=["<3"])

    singles["<3"] = df[tocombine].sum(axis=1)

    final = pd.concat([morethanone, singles], axis=1)

    return final


def general_dt104_mgt_summary(dt104,schemes,outfolder,name,cc_or_st):
    df = pd.DataFrame()
    c=0
    for scheme in schemes:
        # scheme = "id_string_" + scheme

        df2 = dt104[scheme].fillna(0).value_counts().to_frame().T
        # print(df2)
        final = add_singletones(df2).T

        # print(df2)
        if c==0:
            df = final
        else:
            df = pd.merge(df,final,how="outer",left_index=True,right_index=True).fillna(0)
        c+=1
    # df.to_excel(writer, 'overall')

    # print(df)
    plt2 = df.T.plot.bar(stacked=True, figsize=(10, 10), width=0.8,legend=False,fontsize=12).get_figure()
    # plt2.xlabel("Country")
    # plt2.ylabel("Strain count")
    # plt2.rcParams.update({'font.size': 12})
    plt2.tight_layout()

    plt2.savefig(outfolder+name+"_"+cc_or_st+"_counts_nostring.pdf")


hgt = "/Users/michaelpayne/Documents/UNSW/OneDrive - UNSW/HGT_Paper/data/10k_run_initial_29-8-18/MGT_stm4_hgt.txt"
meta = "/Users/michaelpayne/Documents/UNSW/OneDrive - UNSW/HGT_Paper/data/old/10krun_final_21-8-18/10k_full_metadata.xlsx"
outfolder = "/Users/michaelpayne/Documents/UNSW/OneDrive - UNSW/HGT_Paper/data/10k_run_initial_29-8-18/graphs/"

annotated_HGT = ""#"/Users/michaelpayne/Documents/UNSW/OneDrive - UNSW/HGT_Paper/data/10k_run_initial_29-8-18/MGT_stm4_initial.xlsx"#"/Users/michaelpayne/Documents/UNSW/OneDrive - UNSW/HGT_Paper/data/10krun_initial_27-8-18/10k_initial.txt"

schemes = ["MGT"+str(x) for x in range(1,10)]
schemes2 = ["ac_MGT"+str(x) for x in range(2,10)]
print(schemes)

all = data_metadata_import_merge(hgt,meta,annotated_HGT)

all[['Collection Year']] = all[['Collection Year']].replace("nan", '0', regex=True)

all[['Collection Year']] = all[['Collection Year']].astype(int)

for i in range(len(schemes)):
    nscheme = "id_string_" + schemes[i]
    # print(nscheme)
    if i == 0:
        all[nscheme] = all[schemes[i]]
        # print(all[nscheme])
    else:
        all[nscheme] = all["id_string_"+schemes[i-1]] + "-" +all[schemes[i]]
        # print(all["id_string_"+schemes[i-1]])
        # print(all[schemes[i]])
        # print(all["id_string_"+schemes[i-1]] + "-" +all[schemes[i]])
        # sl(0.3)

for i in range(len(schemes2)):
    nscheme = "id_string_" + schemes2[i]
    if i == 0:
        all[nscheme] = all[schemes2[i]]
    else:
        all[nscheme] = all["id_string_"+schemes2[i-1]] + "-" + all[schemes2[i]]
print(all.columns)



# for i in range(len(schemes)):
#     nscheme = "id_string_" + schemes[i]
#     changedict = get_inconsistent(all,nscheme)
#     all[nscheme] = all[nscheme].replace(changedict)



# print(all["id_string_MGT5"])

# print(all.describe())
mather = all.loc[all['substudy']=="mathers"]
print(mather.shape)
leekit = all.loc[all['substudy']=="leekit"]
dt104 = all.loc[all['Initial_experiment']=="DT104"]
nz = all.loc[all['Initial_experiment']=="NZ"]
# writer = ExcelWriter('/Users/michaelpayne/Documents/UNSW/OneDrive - UNSW/HGT_Paper/figures/dt104_combined.xlsx')

leekit_country(leekit,schemes,outfolder)
all_source(mather,schemes,outfolder)
schemes = ["MGT"+str(x) for x in range(1,10)]
mather_UK_time(mather,schemes,outfolder)

# sns.set_palette(sns.color_palette("tab10", 9))
#
major_ST_over_time(dt104,schemes,outfolder,"DT104","st")
major_ST_over_time(nz,schemes,outfolder,"DT160","st")
general_dt104_mgt_summary(dt104,schemes,outfolder,"DT104","st")
general_dt104_mgt_summary(nz,schemes,outfolder,"DT160","st")
# writer.save()
