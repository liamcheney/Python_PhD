import glob
from time import sleep as sl

fail_ls = list(open('/Users/liam/Desktop/failstrains.txt','r').read().splitlines())

for filename in glob.iglob('/Users/liam/Desktop/quast_results/ER*'):
    genome = filename.split('\t')[-1]
    file = open(filename ,'r').read().splitlines()
    for line in file[1:]:
        col = line.split('\t')
        genome = col[0]
        contigs = col[1]
        largest_cont = col[14]
        length = col[15]
        gc = col[16]
        n50 = col[17]
        gene_no = col[22]
        if (genome + ".fasta") in fail_ls:
            print(genome + '\t' + contigs + '\t' + largest_cont + '\t' + length + '\t' + gc + '\t' + n50 + '\t' + gene_no)
            # if 3800000 < int(length) < 4200000:
            #     print(genome, contigs, largest_cont, length, gc, n50, gene_no)