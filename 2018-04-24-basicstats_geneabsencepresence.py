from time import sleep as sl

used to find the average gene number of isolates
infile = open('/Users/liam/Desktop/gene_presence_absence.csv','r').read().splitlines()
839 gnomes total

for i in infile[0:]:
    gcount = 0
    missingcount = 0
    col = i.split(',')
    orth = col[0]
    for j in col:
        if '::' not in j:
            missingcount = missingcount + 1
    gcount = 839 - missingcount
    # print("Gene = " + str(orth) + "     Missing = " + str(missingcount) + "     Carries = " + str(gcount))

misslist = []
col1 = infile[1].split(',')[14:]
print(col1)
if "::" not in col1:
    misslist.append(infile[0])

for k in misslist:
    print(k)
for j in col1:
    if '::' not in j:
        print(j)
        sl(0)

#used to find the presence of absence of an ortholog for a strain
import pandas as pd
desired_width = 500
pd.set_option('display.width', desired_width)

df = pd.read_csv('/Users/liam/Desktop/gene_presence_absence.csv', header=0)

# will determine the amount of strains analysed in the roary file. Below will suggest the amount of strains to join.
straincount = df.shape[1] - 14
top99 = round(straincount - (straincount * 0.01), 0)
#
print(straincount)
print(top99)
#
df['No. isolates'] = pd.to_numeric(df['No. isolates'])
df['Avg sequences per isolate'] = pd.to_numeric(df['Avg sequences per isolate'])
#
df1 = df.loc[df['Avg sequences per isolate'] == 1]

# will iterate over all columns in the csv from 14 onwards, then count the occurence of each gene inside each sub list.
# can tell if a strain is missing certain genes
colcount = 14
runcount = 0
for lab in df.iloc[0,14:]:
    strain = df.columns.values[colcount]
    orthcount = 0
    subdf = df.iloc[ : , colcount]
    colcount = colcount + 1
    for i in subdf.items():
        i = str(i)
        if "ERR" in i:
            orthcount = orthcount + 1
        elif "GCF" in i:
            orthcount = orthcount + 1
        elif "GCA" in i:
            orthcount = orthcount + 1
        elif "cds" in i:
            orthcount = orthcount + 1
    print(str(strain) + "\t" + str(orthcount))
