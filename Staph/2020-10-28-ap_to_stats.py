
import pandas as pd
from time import sleep as sl

def count_values_in_range(series):
    return series.str.match(r'-\d+_\d+')

def per_strain(inp):
    inap = pd.read_csv(inp,sep="\t",dtype=str)
    inap['zeros'] = inap.iloc[:, 3:].isin({'0'}).sum(1)
    inap["negs"] = inap.apply(func=lambda row: count_values_in_range(row), axis=1).sum(1).astype(int)
    outap = inap[["Strain",'zeros','negs']]
    outap.to_csv("/Users/liamcheneyy/Desktop/USA300-NAR_per_strains_alleles.txt",sep="\t",index=False)

def per_locus(inp):
    inap = pd.read_csv(inp,sep="\t",dtype=str)
    locusinfo = inap.iloc[:, 3:].T
    locusinfo['zeros'] = locusinfo.isin({'0'}).sum(1)
    locusinfo["negs"] = locusinfo.apply(func=lambda row: count_values_in_range(row), axis=1).sum(1).astype(int)
    outaploc = locusinfo[['zeros','negs']]
    outaploc.to_csv("/Users/liamcheneyy/Desktop/USA300-NAE_per_locus-alleles.txt",sep="\t")

indata = "/Users/liamcheneyy/Desktop/USA300-NAR-alleles.txt"

per_strain(indata)
per_locus(indata)