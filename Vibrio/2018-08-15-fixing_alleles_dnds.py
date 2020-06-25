from Bio import SeqIO
from collections import defaultdict
from Bio.Seq import translate

inseqs = SeqIO.parse("/Users/liam/Desktop/Vibrio_D3A_alleles_ref.fasta","fasta")
MGT8 = open("/Users/liam/Desktop/species_loci.txt","r").read().splitlines()


c1 = 0
c2 = 0
c3 = 0
cmissed = 0
for i in inseqs:
    s = str(i.seq)
    prot1 = translate(s)
    prot2 = translate(s[1:])
    prot3 = translate(s[2:])
    gene = str(i.id).split(":")[0]
    if gene in MGT8:
        offset = 0
        exist_alleles = SeqIO.parse("/Users/liam/Desktop/Vibrio_D3A_alleles/"+gene+".fasta","fasta")
        if prot1.count("*") == 0 or prot1[-1] == "*":
            c1+=1
            offset = 0
        elif prot2.count("*") == 0 or prot2[-1] == "*":
            c2+=1
            offset = 1
        elif prot3.count("*") == 0 or prot3[-1] == "*":
            c3+=1
            offset = 2
        else:
            offset = -1
            cmissed +=1

        if offset >= 0:
            out = []
            start_list = []
            tga_list = []
            taa_list = []
            tag_list = []
            for al in exist_alleles:
                al.seq = al.seq[offset:]
                sequence = al.seq
                if "N" not in al.seq and "-" not in al.seq:
                    protein = translate(al.seq)

                    start_list.append(protein.count("M"))
                    start_dict = defaultdict(int)
                    for element1 in start_list:
                        start_dict[element1] += 1
                    start_dict = sorted(start_dict.items(), key=lambda x: x[1], reverse=True)
                    start_amount = start_dict[0][0]

                    tga_list.append(sequence.count('TGA'))
                    tga_dict = defaultdict(int)
                    for element2 in tga_list:
                        tga_dict[element2] += 1
                    tga_dict = sorted(tga_dict.items(), key=lambda x: x[1], reverse=True)
                    tga_amount = tga_dict[0][0]

                    taa_list.append(sequence.count('TAA'))
                    taa_dict = defaultdict(int)
                    for element3 in taa_list:
                        taa_dict[element3] += 1
                    taa_dict = sorted(taa_dict.items(), key=lambda x: x[1], reverse=True)
                    taa_amount = taa_dict[0][0]

                    tag_list.append(sequence.count('TAG'))
                    tag_dict = defaultdict(int)
                    for element4 in tag_list:
                        tag_dict[element4] += 1
                    tag_dict = sorted(tag_dict.items(), key=lambda x: x[1], reverse=True)
                    tag_amount = tag_dict[0][0]

                    if (protein.count("M") == start_amount) and \
                            (sequence.count("TGA") == tga_amount and (sequence.count("TAA") == taa_amount) \
                                and (sequence.count("TAG") == tag_amount)):
                        out.append(al)

            if len(out) > 2:
                SeqIO.write(out,"/Users/liam/Desktop/Vibrio_D3A_alleles/codons/"+gene+"_codons.fasta","fasta")
                print(gene,len(out))
            else:
                print("Not enough alleles for loci : " + gene)
