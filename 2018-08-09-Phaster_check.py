MGT_loci = open('/Users/liam/Desktop/Vibrio_D3A_allele_locations.txt','r').read().splitlines()
phast = open('/Users/liam/Desktop/phast.csv').read().splitlines()

# out = open('E:/2018/2018-05-11-Salmonella_enteritidis/test/08-08-07-phast_phage_out.txt','w')
################################################### to get the results of phaster
for a in MGT_loci:
   col = a.split('\t') # col = a.split(',')
   start = int(col[1])
   end = int(col[2])
   locus_tag = col[0]
   for b in phast[1:]:
       col2 = b.split('\t')
       restart = int(col2[2])
       reend = int(col2[3])
       if start <= restart <= end or start <= reend <= end:
           output = col[0], col[1], col[2],col[-1],",".join(col2)
           output = str(output)
           print(output)
