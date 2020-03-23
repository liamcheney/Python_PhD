from time import sleep as sl
import pandas as pd

import sys


def split_hgts(hgtfile,outf):
    hgt = open(hgtfile,"r").read().splitlines()
    if outf:
        hgt_out = open(hgtfile.replace(".txt","_tab.txt"),"w")

    hgt = [x.split("\t") for x in hgt]


    schemes = ['MGT2', 'MGT3', 'MGT4', 'MGT5', 'MGT6', 'MGT7', "MGT8", "MGT9","MGT92","MGT95","MGT910"]

    main_schemes = schemes[:-3]
    extras = schemes[-3:]
    header = []
    head = 0
    rows = []
    strains = []
    for i in hgt:
        if head == 0:
            # print(i)
            head =1
            header = [i[0]]
            for nm in i[1:-3]:
                nm = [nm.replace("mgt",x) for x in main_schemes]
                header += nm
            header += i[-3:]
            if outf:
                hgt_out.write("\t".join(header)+"\n")

        else:
            row = [i[0]]
            strains.append(i[0])
            for nm in i[1:-3]:
                nm = nm.split("-")
                row += nm
            row += i[-3:]
            rows.append(row)
            if outf:
                hgt_out.write("\t".join(row) + "\n")
    df = pd.DataFrame(rows,columns=header,index=strains)
    cols = [c for c in df.columns if "sm" not in c]
    df = df[cols]
    df["MGT91"] = df['ac_MGT9']
    if outf:
        hgt_out.close()

    return df

def main():
    inp = '/Users/liamcheneyy/Desktop/Vibrio_D3A_hgt.txt'
    outptf = '/Users/liamcheneyy/Desktop/Vibrio_tab_D3A_hgt.txt'
    split_hgts(inp,outptf)

if __name__ == "__main__":
    main()