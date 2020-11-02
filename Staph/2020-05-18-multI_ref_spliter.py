import argparse
from time import sleep as sl
from Bio import SeqIO

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args


def main():
    args = parseargs()


    seqs = SeqIO.parse('/Users/liamcheneyy/Desktop/saureus_5/Mgtfi_ref/ref_Staphylococcus_alleles/refs.fasta', 'fasta')
    save_list = []
    for records in seqs:

        locus = records.id.split(':')[0]

        if locus not in save_list:
            save_list.append(locus)

    seqx = SeqIO.to_dict(SeqIO.parse('/Users/liamcheneyy/Desktop/saureus_5/Mgtfi_ref/ref_Staphylococcus_alleles/refs.fasta', 'fasta'))

    save_dict = {}
    for el in save_list:
        save_dict[el] = {}
        for x in seqx:
            if el in x:
                save_dict[el][x] = seqx[x]

    for k,v in save_dict.items():
        with open('/Users/liamcheneyy/Desktop/saureus_5/Mgtfi_ref/ref_Staphylococcus_alleles/' + k + '.fasta','w') as out:
            for i in v:
                out.write('>' + i + '\n')
                out.write(str(v[i].seq) + '\n')

if __name__ == '__main__':
    main()
