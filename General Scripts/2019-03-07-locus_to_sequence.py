from Bio import SeqIO
from time import sleep as sl

infile = open("/Users/liamcheneyy/Desktop/Mgtfi_ref/lociLocationsInRef.txt").read().splitlines()
reference_path = "/Users/liamcheneyy/Desktop/ref/GCA_000006745.1.fasta"
outfile_path = "/Users/liamcheneyy/Desktop/alleles/"

for line in infile:
    col = line.split('\t')
    locus = col[0]
    start = int(col[1]) - 1
    end = int(col[2])
    strand = col[-2]
    chr = col[-1]
    with open(outfile_path + locus + ".fasta", 'w') as out:
        for record in SeqIO.parse(reference_path, "fasta"):
            if chr == record.id:
                seq = record.seq[start:end]

                if strand == '-':
                    seq = record.seq[start:end].reverse_complement()

                out.write(">" + locus + ":1" + '\n')
                out.write(str(seq) + '\n')
                print(locus, strand)


