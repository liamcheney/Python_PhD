from time import sleep as sl
import glob
import pandas as pd
import random
import sys

def too_close(poslists,newpos,limit, chro):
    """

    :param poslists: list of start positions
    :param newpos: start position of gene to test
    :param limit: minimum distance cutoff
    :return: True if within limit, false if further
    """
    poslists1 = [x for x in poslists[chro]]
    poslists1 = map(int,poslists1)
    newpos = int(newpos)
    for i in poslists1:
        if -1*limit <= i-newpos <= limit:
            return True
    return False

input_path = "/Users/liamcheneyy/Desktop/genes_preference_selector.txt"
out_path = input_path.split('/')[1:-1]
out_path = '/' + '/'.join(v for v in out_path)

info = pd.read_csv(input_path, sep="\t").fillna("none")

filt = pd.DataFrame(info)

##assign preferences starting at 13 and overwriting for each more resticted level
filt["pref"] = 11
print("pref 11")
print(filt[filt["pref"]==11]["pref"].count())

filt.loc[(filt["In Species Core"].astype(str).str.contains("T")) & (filt["Seventh Negative Counts"]<=50) & (filt["Seventh Zero Counts"]<=50) & (filt["Species Negative Counts"]<=65) & (filt["Species Zero Counts"]<=65) & (filt["pref"]<=11), "pref"] = 10
print("pref 10")
print(filt[filt["pref"]==10]["pref"].count())

filt.loc[(filt["In dS Ninety Percen"].astype(str).str.contains("T")) & (filt["pref"]<=10), "pref"] = 9
print("pref 9")
print(filt[filt["pref"]==9]["pref"].count())

filt.loc[(filt["In Gene Size Ninety Percen"].astype(str).str.contains("T")) & (filt["In dNdS Ninety Percen"].astype(str).str.contains("T")) & (filt["pref"]<=9), "pref"] = 8
print("pref 8")
print(filt[filt["pref"]==8]["pref"].count())

filt.loc[(filt["In Bio Cycle Excluded"].astype(str).str.contains("F")) & (filt["In Patric Excluded"].astype(str).str.contains("F")) & (filt["pref"]<=8), "pref"] = 7
print("pref 7")
print(filt[filt["pref"]==7]["pref"].count())

filt.loc[(filt["Has tandem repeats"].astype(str).str.contains("F")) & (filt["Has homopolymers"].astype(str).str.contains("F")) & (filt["pref"]<=7), "pref"] = 6
print("pref 6")
print(filt[filt["pref"]==6]["pref"].count())

filt.loc[(filt["Has TMH"].astype(str).str.contains("F")) & (filt["Has Signal Peptides"].astype(str).str.contains("F")) & (filt["Has pred-TAT"].astype(str).str.contains("F")) & (filt["pref"]<=6), "pref"] = 5
print("pref 5")
print(filt[filt["pref"]==5]["pref"].count())

filt.loc[(filt["In Phasta"].astype(str).str.contains("F")) & (filt["In T3SS"].astype(str).str.contains("F")) & (filt["pref"]<=5), "pref"] = 4
print("pref 4")
print(filt[filt["pref"]==4]["pref"].count())

filt.loc[(filt["In dS Fifty Percen"].astype(str).str.contains("T")) & (filt["In dNdS Fifty Percen"].astype(str).str.contains("T")) & (filt["pref"]<=4), "pref"] = 3
print("pref 3")
print(filt[filt["pref"]==3]["pref"].count())

filt.loc[(filt["In Genus Core"].astype(str).str.contains("T")) & (filt["pref"]<=3), "pref"] = 2
print("pref 2")
print(filt[filt["pref"]==2]["pref"].count())

filt.loc[(filt["Seventh Negative Counts"]==0) & (filt["Seventh Zero Counts"]==0) & (filt["Species Negative Counts"]==0) & (filt["Species Zero Counts"]==0) & (filt["pref"]<=2), "pref"] = 1
print("pref 1")
print(filt[filt["pref"]==1]["pref"].count())


#choosing genes in MGT2 and MGT3
# filt["pref"] = 2
# print("pref 2")
# print(filt[filt["pref"]==2]["pref"].count())
#
# filt.loc[(filt["Seventh Negative Counts"]<=50) & (filt["Seventh Zero Counts"]<=50) & (filt["Species Negative Counts"]<=60) & (filt["Species Zero Counts"]<=60) & (filt["In Spc All Changing"].astype(str).str.contains("T")) & (filt["pref"]<=2), "pref"] = 1
# # filt.loc[(filt["Number of Negative Counts"]<=0) & (filt["Number of Zero Counts"]<=0) & (filt["In All Spc Alle"].astype(str).str.contains("T")) & (filt["In Species Core"].astype(str).str.contains("T")) & (filt["pref"]<=2), "pref"] = 1
# print("pref 1")
# print(filt[filt["pref"]==1]["pref"].count())

# sys.exit()

filt.to_csv('/Users/liamcheneyy/Desktop/filt_all_genes_hgt.csv', index=False)

#for MGT2 and MGT3
random.seed(561618)

#scheme target sizes
# target_sizes = {'MGT2':10329,'MGT3':51644, 'MGT4':103287}
#
# ##scheme lowest allowed loci preference numbers
# preflimit = {'MGT2':1,'MGT3':1, 'MGT4':1}
#
# ##scheme smallest distance allowed between loci
# distlimit = {'MGT2':20000,'MGT3':10000, 'MGT4':500}

#For MGT4 onwards
target_sizes = {'MGT5':206575,'MGT6':516437,'MGT7':1032875}

#scheme lowest allowed loci preference numbers
preflimit = {'MGT5':5,'MGT6':8,'MGT7':10}

#scheme smallest distance allowed between loci
distlimit = {'MGT5':3500,'MGT6':1000,'MGT7':0}

#initialise outputs dict

#read in manually assigned genes
assigned_genes = {}
for filename in glob.iglob('/Users/liamcheneyy/Desktop/manually_assigned_genes/*'):
    mgt = filename.split('/')[-1].strip('.txt')
    file = open(filename,'r').read().splitlines()
    assigned_genes[mgt] = file

#creating output
outputs = assigned_genes


#list of used genes
donegenes = []
for key, value in assigned_genes.items():
    for i in value:
        donegenes.append(i)

#length of total input genes
input_length = {}
for key, value in assigned_genes.items():
    length = 0
    for el in value:
        gene_length = int(filt[filt['Locus Tag'] == el]['Length'].values)
        length = length + gene_length
    input_length[key] = length

#list of gene start pos
startposs = {2:[], 1:[]}
for key, value in assigned_genes.items():
    for el in value:
        if 'VCA' not in el:
            start = int(filt[filt['Locus Tag'] == el]['Start'].values)
            startposs[1].append(start)
        if el[0:3] == 'VCA':
            start = int(filt[filt['Locus Tag'] == el]['Start'].values)
            startposs[2].append(start)

#dict of genes to preferences
prefassigns = {}
for el in donegenes:
    prefassigns[el] = 1

#dict of distance to list of loci that are less that that distance
toclose_genes = {x:[] for x in range(0,20001,1000)}
# print(toclose_genes)

#preference number
prefno = 1
#length limit
limno = 20000

#get initial list of loci names from preference 1
pref_loci = list(filt[filt["pref"]==prefno]["Locus Tag"].get_values())

#make a dataframe of oly preference 1
pref_df = pd.DataFrame(filt[filt["pref"]==prefno])
# print(pref_loci)
#for each MGT scheme
for i in sorted(target_sizes.keys()):
    print(i)
    outputs[i] = assigned_genes[i]
    totlen = input_length[i]
    locino = len(assigned_genes[i])
    # outputs[i] = []
    # totlen = 0
    # locino = 0

    #while current running total of scheme genes is less than the target size
    while totlen < target_sizes[i]:
        #if there are any genes left in the pref_loci list
        if len(pref_loci) > 0:
            #get number of genes
            geneno = len(pref_loci)
            #get a random number from 0 to above number
            a = random.randint(0, geneno-1)

            #get gene name
            ids = pref_loci[a]
            #get gene info from df
            gene = filt.loc[filt["Locus Tag"]==ids]
            pos = int(gene["Start"])
            length = int(gene["Length"])
            chro = gene["Chromosome"].values[0]
            # if chro not in startposs:
            #     startposs[chro] = []
            #check if new gene is too close to existing picks given current limit (limno) using start position of new gene
            to_close_filt = too_close(startposs,pos,limno, chro)
            ##if the new gene is not already picked and is not already noted as too close fr current limit
            if ids not in donegenes and ids not in toclose_genes[limno]:
                #if the new gene is not too close
                if to_close_filt == False:
                    #add gene to done list, outputs for current MGT, start position list, preference recording, total length so far for MGT, number of loci, remove gene from current preference list
                    donegenes.append(ids)
                    outputs[i].append(ids)
                    startposs[chro].append(pos)
                    prefassigns[ids] = prefno
                    totlen += length
                    locino+=1
                    pref_loci.remove(ids)
                else:
                    ## if gene is too close add to to_close_genes dictionary for current limit
                    toclose_genes[limno].append(ids)
                    ## remove gene from current preference list
                    pref_loci.remove(ids)
            else:
                #if gene is already in done or toclose remove from current preference list
                pref_loci.remove(ids)
        #if there are no more genes in current preference and the preference number is under the limit for the MGT
        elif prefno < preflimit[i]:
            # print("HELLLOOOO2")
            #increase the preference no by 1
            prefno +=1
            # get loci for this preference
            pref_loci = list(filt[filt["pref"]==prefno]["Locus Tag"].get_values())
            pref_df = pd.DataFrame(filt[filt["pref"] == prefno])
            # print("preference_no: "+str(prefno))
        # if length limit is  is more than cutoff for min distance
        elif limno > distlimit[i]:
            # reset preference to 1 and reset pref_loci to preference 1
            prefno = 1
            pref_loci = list(filt[filt["pref"]==prefno]["Locus Tag"].get_values())
            #reduce distance limit to 1000 less
            limno -= 1000
            #if limit is more than 0
            if limno > 0:
                #add lower cuttoff genes to current - because genes in 6000 cutoff chould be included in 7000 cutoff etc
                toclose_genes[limno] += toclose_genes[limno-1000]

            # print(toclose_genes)
            print("distance limit: "+str(limno),len(toclose_genes[limno]))
        else:
            break
    #test for getting loci when not enough sequence is included
    if totlen < target_sizes[i]:
        print("NOT ENOUGH SEQ!!\nNeed: "+str(target_sizes[i]))
        print(i,totlen,locino)

        sl(100000)
    #reset distance limit for each scheme
    limno = 20000


# output writing
# print(len(toclose_genes))
# print(len(donegenes))
outfolder = '/Users/liamcheneyy/Desktop/loci/'
outsummary = open(outfolder + "/all_schemes_loci.txt","w")
out_alleles = '/Users/liamcheneyy/Desktop/refonly_allelic_profiles/'
for i in outputs:
    outf = open(outfolder +"/"+i+"_gene_accessions.txt","w")
    print(i,len(outputs[i]))
    for gene in outputs[i]:
        outsummary.write("{}\t{}\t{}\n".format(i,gene,prefassigns[gene]))
    outf.write("\n".join(outputs[i]))
    outf.close()
outsummary.close()

#writeout allelic profiles
for i in outputs:
    with open(out_alleles +"/"+i+"_gene_profiles.txt","w") as out:
        out.write('ST' + '\t' + 'dST' + '\t')
        for gene in outputs[i]:
            out.write(gene + '\t')
        out.write('\n')
        out.write('1' + '\t' + '0' + '\t')
        for gene in outputs[i]:
            out.write('1' + '\t')

