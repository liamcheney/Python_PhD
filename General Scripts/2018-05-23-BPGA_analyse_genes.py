from time import sleep as sl
from Bio import SeqIO

##search the sequences from BPGA core seq file output in references sequence file. Convert to protein ID.
org_dict = {}
org1_seq = SeqIO.parse('/Users/liam/Desktop/analysis_D1a_BPGA_45%true/Org1_ref_core_seq.txt', 'fasta')
for i in org1_seq:
    org_dict[str(i.id)] = str(i.seq)

ref_dict = {}
ref_seq = SeqIO.parse('/Users/liam/Desktop/analysis_D1a_BPGA_45%true/GCA_000006745.1.fas', 'fasta')
for j in ref_seq:
    ref_dict[str(j.seq)] = str(j.id)

bpga_cor_pro_id = []
for k in org_dict.values():
    bpga_cor_pro_id.append(ref_dict[k])

## reading in reference .gpff file to convert the protein IDs into VC locus tags. Also find cds region and search in ref .gff to get 3 columns of protein, locus and cds tags.
from time import sleep as sl
from Bio import SeqIO

infile = open('/Users/liam/Desktop/analysis_D1a_BPGA_45%true/GCA_000006745.1_ASM674v1_protein.gpff', 'r').readlines()
outfile = open('/Users/liam/Desktop/analysis_D1a_BPGA_45%true/protein_locus_cds.txt', 'w')

my_dict = {}
for i in infile:
    if "VERSION" in i:
        protein_ID = i.split('     ')[1]
        protein_ID = protein_ID.rstrip('\n')
        my_dict[protein_ID] = 'x'
    if "locus_tag" in i:
        locus_tag = i.split('                     /locus_tag="')[1].rstrip('"\n')
        ocus_tag = locus_tag.replace('VC_A', 'VCA').replace('VC_', 'VC')
        my_dict[protein_ID] = locus_tag
    # if '/coded_by=' in i:
    #     cds_cor = i.split(':')[1].rstrip(')"\n').split('..')
    #     ref_gff_in = open('/Users/liam/Desktop/analysis_D1a_BPGA_45%true/GCF_000006745.1_ASM674v1_genomic.gff','r').readlines()
        # for line in ref_gff_in:
            # if str(cds_cor[1]) in line:

#searching the protein_ID from bgpga in the reference strain and converting to locus_tags.
# for l in bpga_cor_pro_id:
#     print(my_dict[l])

#taking cds to locus info from excel file
infile1 = open('/Users/liam/Desktop/ref_cds_to_locus.txt','r').read().splitlines()

ref_info_dict = {}
for i in infile1:
    col = i.split('\t')
    ref_info_dict[col[0]] = col[1]

#read in roary cds and conert to locus_tag using ref_info_dict
roary_cds_in = open('/Users/liam/Desktop/roary_cds.txt','r').read().splitlines()

for cds in roary_cds_in:
    print(ref_info_dict[cds])