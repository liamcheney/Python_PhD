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

    # infiles = open("/Users/liamcheney/Desktop/all_ncbi_outputs.csv").read().split('\n\n')
    #
    # save_dict = {}
    # save_list = []
    # for file in infiles:
    #     infile = file.split('\n')
    #     for line in infile[1:]:
    #         col = line.split(',')
    #         gene = col[4]
    #         if gene not in save_list:
    #             save_list.append(gene)

    # print(len(save_list))
    # for i in save_list:
    #     print(i)

    # save_dict = {}
    # for filename in glob.iglob('/Users/liamcheney/Desktop/outs/*_out.txt'):
    #     infile = open(filename).read().splitlines()
    #     acc = filename.split('/')[-1].split('_')[0]
    #     intact = 0
    #     partial = 0
    #     called = 0
    #     print(acc)
    #
    #     for line in infile:
    #
    #         if "intact:" in line:
    #             intact = int(line.split(':')[-1].strip())
    #         if "partial:" in line:
    #             partial = int(line.split(':')[-1].strip())
    #         if "called:" in line:
    #             called = int(line.split(':')[-1].strip())
    #
    #     total_pass = intact + partial + called
    #     failed =  1860 - total_pass
    #     save_dict[acc] = {'total':total_pass, 'fail':failed}
    #
    # with open('/Users/liamcheney/Desktop/all_out11.csv','w') as out:
    #     for k,v in save_dict.items():
    #         out.write(k + ',')
    #         for x in v:
    #             out.write(str(save_dict[k][x]) + ',')
    #         out.write('\n')


    # ##check how gene went across dataset
    # infile = open('//Users/liamcheneyy/Desktop/all_alleles.txt').read().splitlines()
    # # print(len(infile))
    #
    # #make list of genes
    # save_dict = {}
    # for line in infile[0:1]:
    #     col = line.split(' ')
    #    # acc=col[0].split('/')[-1].split('_')[0]
    #     for j in col[4:]:
    #         gene = j.split(':')[0].strip('>')
    #         save_dict[gene] = {'pass':0, 'fail':0}
    #
    # #iterate over
    # for line in infile:
    #     col = line.split(' ')
    #
    #     acc=col[0].split('/')[-1].split('_')[0]
    #
    #     for j in col[4:]:
    #         gene = j.split(':')[0].strip('>')
    #         if ':0' in j:
    #             save_dict[gene]['fail'] = save_dict[gene]['fail'] + 1
    #         # else:
    #         #     save_dict[gene]['pass'] = save_dict[gene]['pass'] + 1
    # #
    # for k,v in save_dict.items():
    #     print(k,v['pass'],v['fail'], sep='\t')
    # #
    # # ##check per strain zeros
    # # save_dict = {}
    # # for line in infile:
    # #     col = line.split(' ')
    # #
    # #     acc=col[0].split('/')[-1].split('_')[0]
    # #
    # #     print(acc, line.count(':0'))

    # ##count negs per gene
    # path = "/Users/liamcheneyy/Desktop/all/"
    #
    # save_dict = {}
    # for file in glob.iglob(path + "/*alleles.fasta"):
    #     seqs = SeqIO.parse(file, 'fasta')
    #     name=file.split('/')[-1].split('.')[0]
    #     for records in seqs:
    #
    #         locus = records.id.split(':')[0]
    #
    #         if locus not in save_dict.keys():
    #             save_dict[locus] = 0
    #
    #         if (':1' in records.id) or (':2' in records.id) or (':3' in records.id) or (':4' in records.id):
    #             save_dict[locus] = save_dict[locus] + 1
    #
    # for k,v in save_dict.items():
    #     print(k, v, sep='\t')

# infile = open('/Users/liamcheneyy/Desktop/MLST_test.txt').read().splitlines()
    #
    # with open('/Users/liamcheneyy/Desktop/MLST_pass.txt','w') as out:
    #     for line in infile:
    #         col = line.split('\t')
    #         if col[1] != '-':
    #             if col[2] == '-':
    #                 if ('(-)' not in line) and ('?' not in line):
    #                     print(line)

    # #FIX OUTPUTS TO REMOVE POOR ASSEMBLES
    # infiles = open('/Users/liamcheneyy/Desktop/all_resfinder.csv').read().split('./')
    # remove = open('/Users/liamcheneyy/Desktop/remove.txt').read().splitlines()
    #
    # with open('/Users/liamcheneyy/Desktop/90_all_alleles_out_fix.csv','w') as out:
    #
    #     for file in infiles[1:-1]:
    #         infile = file.splitlines()
    #         acc = infile[0].split('/')[-1].split('_')[0]
    #         if acc not in remove:
    #             out.write('./' + file)


    # with open('/Users/liamcheneyy/Desktop/90_all_alleles_out_fix.csv','w') as out:
    #     for line in infiles:
    #         col = line.split(',')
    #         acc=col[0].split('_')[0]
    #         print(acc)
    #         # if acc not in remove:
    #         #     out.write(line + '\n')



    # seqs = SeqIO.parse('/Users/liamcheneyy/Desktop/refs.fasta', 'fasta')
    # save_list = []
    # for records in seqs:
    #
    #     locus = records.id.split(':')[0]
    #
    #     if locus not in save_list:
    #         save_list.append(locus)
    #
    # seqx = SeqIO.to_dict(SeqIO.parse('/Users/liamcheneyy/Desktop/refs.fasta', 'fasta'))
    #
    # save_dict = {}
    # for el in save_list:
    #     save_dict[el] = {}
    #     for x in seqx:
    #         if el in x:
    #             save_dict[el][x] = seqx[x]
    #
    # for k,v in save_dict.items():
    #     with open('/Users/liamcheneyy/Desktop/outs/' + k + '.fasta' ,'w') as out:
    #         for i in v:
    #             out.write('>' + i + '\n')
    #             out.write(str(v[i].seq) + '\n')

    ###go over allele profiles and get summary stats
    genome_zero_neg_dict = {}
    gene_zero_dict = {}
    for filename in glob.iglob('/Users/liamcheneyy/Desktop/alleles/*fasta'):
        bios = SeqIO.parse(filename,'fasta')
        accession = filename.split('/')[-1].split('_')[0]
        print(accession)

        genome_zero_neg_dict[accession] = {'zero':0,'neg':0}

        ##count zeros per genome
        for record in bios:
            gene = record.id.split(':')[0]

            if gene not in gene_zero_dict.keys():
                gene_zero_dict[gene] = {'zero':0, 'intact':0, 'new':0, 'neg':0}

            if ":0" in record.id:
                gene_zero_dict[gene]['zero'] = gene_zero_dict[gene]['zero'] + 1
                genome_zero_neg_dict[accession]['zero'] = genome_zero_neg_dict[accession]['zero'] + 1

            if ":1" in record.id or ":2" in record.id or ":3" in record.id or ":4" in record.id or ":5" in record.id or ":6" in record.id or ":7" in record.id:
                gene_zero_dict[gene]['intact'] = gene_zero_dict[gene]['intact'] + 1

            if ":new" in record.id:
                gene_zero_dict[gene]['new'] = gene_zero_dict[gene]['new'] + 1

            if 'N' in record.seq:
                gene_zero_dict[gene]['neg'] = gene_zero_dict[gene]['neg'] + 1
                genome_zero_neg_dict[accession]['neg'] = genome_zero_neg_dict[accession]['neg'] + 1


    with open('/Users/liamcheneyy/Desktop/alleles/genome_zeros.txt', 'w') as out:
        out.write('Accessions' + '\t' + 'Zero' + '\t' + 'Neg' + '\n')
        for key, value in genome_zero_neg_dict.items():
            out.write(key + '\t')
            for i in value:
                out.write(str(value[i]) + '\t')
            out.write('\n')

    with open('/Users/liamcheneyy/Desktop/alleles/gene_zeros.txt', 'w') as out:
        out.write('Gene' + '\t' + 'Zero' + '\t' + 'Intact' + '\t' + 'New' + '\t' + 'Neg' + '\n')
        for key, value in gene_zero_dict.items():
            out.write(key + '\t')
            for i in value:
                out.write(str(value[i]) + '\t')
            out.write('\n')




        # genome_zero_dict[accession] = genome_zero_conut

if __name__ == '__main__':
    main()