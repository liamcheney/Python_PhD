from time import sleep as sl
import glob
import pandas as pd
import random
import sys

##Adjsutable
def preferences(info):
    filt = pd.DataFrame(info)

    p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11 = 1,2,3,4,5,6,7,8,9,10,11

    ##assign preferences
    filt["pref"] = p11
    print(f"pref p11", filt[filt["pref"] == p11]["pref"].count())

    filt.loc[(filt["alleles_96_percen"].astype(str).str.contains("T")) & (filt["pref"] <= p10+1), "pref"] = p10
    print(f"pref {p1}", filt[filt["pref"] == p10]["pref"].count())

    filt.loc[(filt["dNdS_90_percentile"].astype(str).str.contains("T")) & (filt["length_90_percentile"].astype(str).str.contains("T")) & (filt["Recombination Events"] <= 30) & (filt["pref"] <= p9+1), "pref"] = p9
    print(f"pref {p9}", filt[filt["pref"] == p9]["pref"].count())

    filt.loc[(filt["Homopolymer"].astype(str).str.contains("F")) & (filt["In Phage Region"].astype(str).str.contains("F")) & (filt["Tandem Repeats"].astype(str).str.contains("F")) & (filt["pref"] <= p8+1), "pref"] = p8
    print(f"pref {p8}", filt[filt["pref"] == p8]["pref"].count())

    filt.loc[(filt["alleles_70_percen"].astype(str).str.contains("T")) & (filt["pref"] <= p7+1), "pref"] = p7
    print(f"pref {p7}", filt[filt["pref"] == p7]["pref"].count())

    filt.loc[(filt["Recombination Events"] <= 20) & (filt["Zero_percen"] <= 0.5) & (filt["Negative_percen"] <= 1) & (filt["pref"] <= p6+1), "pref"] = p6
    print(f"pref {p6}", filt[filt["pref"] == p6]["pref"].count())

    filt.loc[(filt["alleles_60_percen"].astype(str).str.contains("T")) & (filt["pref"] <= p5+1), "pref"] = p5
    print(f"pref {p5}", filt[filt["pref"] == p5]["pref"].count())

    filt.loc[(filt["alleles_50_percen"].astype(str).str.contains("T")) & (filt["dNdS_50_percentile"].astype(str).str.contains("T")) & (filt["pref"] <= p4+1), "pref"] = p4
    print(f"pref {p4}", filt[filt["pref"] == p4]["pref"].count())

    filt.loc[(filt["Recombination Events"] <= 10) & (filt["Zero_percen"] <= 0.05) & (filt["Negative_percen"] <= 0.05) & (filt["pref"] <= p3+1), "pref"] = p3
    print(f"pref {p3}", filt[filt["pref"] == p3]["pref"].count())

    filt.loc[(filt["alleles_20_percen"].astype(str).str.contains("T")) & (filt["pref"] <= p2+1), "pref"] = p2
    print(f"pref {p2}", filt[filt["pref"] == p2]["pref"].count())

    filt.loc[(filt["alleles_10_percen"].astype(str).str.contains("T")) & (filt["pref"] <= p1+1), "pref"] = p1
    print(f"pref {p1}", filt[filt["pref"] == p1]["pref"].count())

    #filt.to_csv('/Users/liamcheneyy/Desktop/filt_all_genes_hgt.csv', index=False)

    return filt

def wanted_limits():
    # for MGT2 and MGT3
    random.seed(5432213243)

    ##scheme target sizes
    target_sizes = {'MGT2':17559,'MGT3':35118,'MGT4':70236,'MGT5':175589,'MGT6':351178,'MGT7':702356}

    # scheme lowest allowed loci preference numbers
    preflimit = {'MGT2': 1, 'MGT3': 3,'MGT4': 4, 'MGT5': 5,'MGT6': 9, 'MGT7': 10}

    # scheme smallest distance allowed between loci
    distlimit = {'MGT2': 20000, 'MGT3': 10000,'MGT4': 5000, 'MGT5': 1000,'MGT6': 0, 'MGT7': 0}

    return target_sizes, preflimit, distlimit

##workings
def too_close(poslists,newpos,limit):
    """

    :param poslists: list of start positions
    :param newpos: start position of gene to test
    :param limit: minimum distance cutoff
    :return: True if within limit, false if further
    """

    # print(poslists, newpos, limit)
    poslists1 = [x for x in poslists]
    poslists1 = map(int,poslists1)

    newpos = int(newpos)
    for i in poslists1:
        if -1*limit <= i-newpos <= limit:
            return True
    return False

def initialise_variables(filt):

    # initialise outputs dict
    # creating output
    outputs = {}

    # list of used genes
    donegenes = []

    # list of gene start pos
    startposs = []

    # dict of genes to preferences
    prefassigns = {}

    # dict of distance to list of loci that are less that that distance
    toclose_genes = {x: [] for x in range(0, 80001, 1000)}
    # print(toclose_genes)

    # preference number
    prefno = 1
    # length limit
    limno = 80000

    # get initial list of loci names from preference 1
    pref_loci = list(filt[filt["pref"] == prefno]["Locus Tag"])

    # make a dataframe of noly preference 1
    pref_df = pd.DataFrame(filt[filt["pref"] == prefno])

    return outputs, donegenes, startposs, prefassigns, toclose_genes, prefno, limno, pref_loci, pref_df

def choose_genes(outputs, donegenes, startposs, prefassigns, toclose_genes, prefno, limno, pref_loci, pref_df, target_sizes, preflimit, distlimit, filt):
    # print(pref_loci)
    # for each MGT scheme
    for i in sorted(target_sizes.keys()):
        print("Starting: ", i)
        outputs[i] = []
        totlen = 0
        locino = 0

        # while current running total of scheme genes is less than the target size
        while totlen < target_sizes[i]:
            # if there are any genes left in the pref_loci list
            if len(pref_loci) > 0:
                # get number of genes
                geneno = len(pref_loci)
                # get a random number from 0 to above number
                a = random.randint(0, geneno - 1)

                # get gene name
                ids = pref_loci[a]
                # get gene info from df
                gene = filt.loc[filt["Locus Tag"] == ids]
                pos = int(gene["Start"])
                length = int(gene["Length"])
                # check if new gene is too close to existing picks given current limit (limno) using start position of new gene
                to_close_filt = too_close(startposs, pos, limno)
                # print(ids, pos, to_close_filt)
                # print(startposs)
                ##if the new gene is not already picked and is not already noted as too close fr current limit
                if ids not in donegenes and ids not in toclose_genes[limno]:
                    # if the new gene is not too close
                    if to_close_filt == False:
                        # add gene to done list, outputs for current MGT, start position list, preference recording, total length so far for MGT, number of loci, remove gene from current preference list
                        donegenes.append(ids)
                        outputs[i].append(ids)
                        startposs.append(pos)
                        prefassigns[ids] = prefno
                        totlen += length
                        locino += 1
                        pref_loci.remove(ids)
                    else:
                        ## if gene is too close add to to_close_genes dictionary for current limit
                        toclose_genes[limno].append(ids)
                        ## remove gene from current preference list
                        pref_loci.remove(ids)
                else:
                    # if gene is already in done or toclose remove from current preference list
                    pref_loci.remove(ids)
            # if there are no more genes in current preference and the preference number is under the limit for the MGT
            elif prefno < preflimit[i]:
                # increase the preference no by 1
                prefno += 1
                # get loci for this preference
                pref_loci = list(filt[filt["pref"] == prefno]["Locus Tag"])
                pref_df = pd.DataFrame(filt[filt["pref"] == prefno])
                # print("preference_no: "+str(prefno))
            # if length limit is  is more than cutoff for min distance
            elif limno > distlimit[i]:
                # print(distlimit, i, limno)
                # reset preference to 1 and reset pref_loci to preference 1
                prefno = 1
                pref_loci = list(filt[filt["pref"] == prefno]["Locus Tag"])
                # reduce distance limit to 1000 less
                limno -= 1000
                # if limit is more than 0
                if limno > 0:
                    # add lower cuttoff genes to current - because genes in 6000 cutoff chould be included in 7000 cutoff etc
                    toclose_genes[limno] += toclose_genes[limno - 1000]

                # print(toclose_genes)
                # print("distance limit: " + str(limno), len(toclose_genes[limno]))
            else:
                break


        # test for getting loci when not enough sequence is included
        if totlen < target_sizes[i]:
            print("NOT ENOUGH SEQ!!\nNeed: " + str(target_sizes[i]))
            print(i, totlen, locino)
            sys.exit()

        # reset distance limit for each scheme
        limno = 20000

    return toclose_genes, outputs, prefassigns

def output(toclose_genes, outputs, prefassigns):
    # output writing
    # print(len(toclose_genes))
    # print(len(donegenes))
    outfolder = '/Users/liamcheney/Desktop/loci/'
    outsummary = open(outfolder + "/all_schemes_loci.txt", "w")
    out_alleles = '/Users/liamcheney/Desktop/refonly_allelic_profiles/'
    for i in outputs:
        outf = open(outfolder + "/" + i + "_gene_accessions.txt", "w")
        print(i, len(outputs[i]))
        for gene in outputs[i]:
            outsummary.write("{}\t{}\t{}\n".format(i, gene, prefassigns[gene]))
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

def main():

    input_path = "/Users/liamcheney/Desktop/Preferences.xlsx"

    info = pd.read_excel(open(input_path, 'rb'), sheet_name='Preferences').fillna("none")

    filt = preferences(info)

    target_sizes, preflimit, distlimit = wanted_limits()

    outputs, donegenes, startposs, prefassigns, toclose_genes, prefno, limno, pref_loci, pref_df = initialise_variables(filt)

    toclose_genes, outputs, prefassigns = choose_genes(outputs, donegenes, startposs, prefassigns, toclose_genes, prefno, limno, pref_loci, pref_df, target_sizes, preflimit, distlimit, filt)

    output(toclose_genes, outputs, prefassigns)

if __name__ == '__main__':
    main()
