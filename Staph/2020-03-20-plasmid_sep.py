import argparse
from time import sleep as sl
from Bio import SeqIO
import glob

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args


def main():
    args = parseargs()

    infolder = "/Users/liamcheneyy/Downloads/NCBI_complete_assemblies/GCA*"
    outfolder = "/Users/liamcheneyy/Downloads/NCBI_complete_assemblies/edits/"

    save_dict = {}
    for filename in glob.iglob(infolder):
        accession = filename.split('/')[-1].split('.')[0]
        chromsome = []
        plasmid = []

        for record in SeqIO.parse(filename,"fasta"):
            description = record.description

            if "plasmid" not in description:
                chromsome.append(record)

            elif "plasmid" in description:
                plasmid.append(record)

            else:
                print(accession)

        SeqIO.write(chromsome,outfolder + "/chromosome/" + accession + ".fna","fasta")

        SeqIO.write(plasmid,outfolder + "/plasmids/" + accession + "plasmids.fna","fasta")

if __name__ == '__main__':
    main()
