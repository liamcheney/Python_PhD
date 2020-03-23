from time import sleep as sl
import pprint
import sys
instrain_ls = []
inprof = open('/Users/liamcheneyy/Desktop/MGT9_gene_profiles.txt',"r").read().splitlines()
inassignment = open('/Users/liamcheneyy/Desktop/MGT9_gene_annotation.txt',"r").read().splitlines()
outprof = open('/Users/liamcheneyy/Desktop/all_data_grapetree.txt',"w")
instrain_ls = open('/Users/liamcheneyy/Desktop/strains.txt',"r").read().splitlines()

h=0
st_to_profile = {}


for i in inprof:
    cols = i.split("\t")
    if h==0:
        outprof.write("#Strain\t" + "\t".join([cols[0]]+cols[2:]) + "\n")
        h=1
    else:
        st = ".".join(cols[:2])
        prf = []
        cols = i.split("\t")[2:]
        # outprof.write("\t".join(cols[0]))
        for i in cols:
            allele = i
            if allele == "0":
                prf.append("-")
                # outprof.write("\t-")
            elif "-" in allele:
                n = allele.split("_")[0]
                # outprof.write("\t"+n)
                prf.append(n[1:])
            else:
                prf.append(allele)
                # outprof.write("\t" + allele)
        st_to_profile[st] = prf
pprint.pprint(st_to_profile.keys())
#
# newst = {}
# stoptions = {}
# for fst in st_to_profile:
#     st = fst.split(".")[0]
#     dst = fst.split(".")[1]
#     if dst == 0:
#         newst[st] = st_to_profile[fst]
#     else:
#         if st not in stoptions:
#             stoptions[st] = [fst]
#         else:
#             stoptions[st].append(fst)
# for fst in st_to_profile:
#     st = fst.split(".")[0]
#     dst = fst.split(".")[1]
#     if st not in newst:
#         newst[st] = st_to_profile[stoptions[st][0]]


for i in inassignment[1:]:
    col = i.split("\t")

    strain = col[0]
    if instrain_ls != []:
        if strain in instrain_ls:
            st = col[1]#.split(".")[0]
            if "." not in st:
                st = st+".0"
            if st[0] != '0':
                prof = st_to_profile[st]
                shortst = st.split(".")[0]
                outprof.write("{}\t{}\t{}\n".format(strain,shortst,"\t".join(prof)))
    else:
        st = col[1]  # .split(".")[0]
        if "." not in st:
            st = st + ".0"
        if st != "0.0":
            prof = st_to_profile[st]
            shortst = st.split(".")[0]
            outprof.write("{}\t{}\t{}\n".format(strain, shortst, "\t".join(prof)))


outprof.close()

    # outprof.write("\n")
# outprof.close()
