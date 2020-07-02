import argparse
from time import sleep as sl

def main():
    recom_snps = open('/Users/liamcheneyy/Desktop/norecom.rec').read().splitlines()

    infile = open('/Users/liamcheneyy/Desktop/lociLocationsInRef.txt').read().splitlines()

    save = {}
    for locuses in infile:
        complete_overlaps = 0
        front_partial_overlaps = 0
        end_partial_overlaps = 0
        col = locuses.split('\t')
        locus = col[0]
        locus_s = int(col[1])
        locus_e = int(col[2])
        for recom_regions in recom_snps:
            col1 = recom_regions.split('\t')
            rec_s = int(col1[1])
            rec_e = int(col1[2])
            if (rec_s <= locus_s) and (rec_e >= locus_e):
                complete_overlaps += 1

            elif ((rec_s <= locus_s) and (rec_e <= locus_e and rec_e >= locus_s)):
                front_partial_overlaps += 1

            elif ((rec_e >= locus_e) and (rec_s >= locus_s and rec_s <= locus_e)):
                end_partial_overlaps += 1

        save[locus] = {complete_overlaps, front_partial_overlaps, end_partial_overlaps}
        total_overlaps = complete_overlaps + front_partial_overlaps + end_partial_overlaps
        print(locus, complete_overlaps, front_partial_overlaps, end_partial_overlaps, total_overlaps,sep='\t')




if __name__ == '__main__':
    main()
