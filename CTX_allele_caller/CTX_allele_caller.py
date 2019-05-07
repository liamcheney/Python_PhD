import argparse
import  glob
from os import path
import sys, os
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import multiprocessing
from time import sleep as sl

def input_allele_lengths(set_wd, args):

    gene_size_dict = {}
    for line in open(set_wd + "/allele_lengths.txt",'r').read().splitlines():
        col = line.split('\t')
        gene_size_dict[col[0]] = int(col[1])
    return gene_size_dict

def choose_database(set_wd, args):
    #get list of all databases
    total_db_list = [x.split('/')[-1] for x in glob.glob(set_wd + "/blast_db/*")]

    #will create a list of databases to analyse genomes against, fist existing list then read user input
    using_db_dict = {}
    for i in args.database:

        if i == '4':
            for n in total_db_list:
                using_db_dict[n] = {}

        else:
            nums = i.split(',')
            for j in nums:
                using_db_dict[total_db_list[int(j)]] = {}

    return using_db_dict

def iterate_databases(databases, set_wd, args):

    gene_size_dict = input_allele_lengths(set_wd, args)

    for datab in databases:
        print("Analaying strains use database: " + datab)
        current_blast_db_path = set_wd + "/blast_db/" + datab + "/"
        results_dict = run_genomes(current_blast_db_path, gene_size_dict, args)
        databases[datab] = results_dict

    return databases

def run_genomes(current_blast_db_path, gene_size_dict, args):
    print("Processing genomes")
    results_dict = {}
    for strain in glob.iglob(args.strains_directory + '/*.f*'):
        accession = strain.split('/')[-1].split('.')[0]
        print(accession)
        results_dict[accession] = blast_input_against_db(strain, current_blast_db_path, gene_size_dict, args)

    return results_dict

def blast_input_against_db(strain, current_blast_db_path, gene_size_dict, args):

    cpus = multiprocessing.cpu_count()
    blast_db_path = current_blast_db_path + current_blast_db_path.split('/')[-2] + "_alleles_db"

    #only want top hit
    if args.blast_return == 1:
        blast_string = NcbiblastnCommandline(query=strain,db=blast_db_path, outfmt=6, perc_identity=50, num_threads=cpus, max_target_seqs=1)
        out, err = blast_string()

    #want all hits ##WARNING: many with low lengths
    elif args.blast_return == 2:
        blast_string = NcbiblastnCommandline(task="blastn", query=strain, db=blast_db_path, outfmt=6, perc_identity=50, num_threads=cpus)
        out, err = blast_string()

    #if not blast alignment found
    if out == '' and '1' not in args.database:
        blast_result_dict = {1:'No alignment found. Modify "-r 2" to see low quality hits.'}
        return blast_result_dict

    elif out == '' and '1' in args.database:
        blast_result_dict = {1: 'non-O1/non-O139'}
        return blast_result_dict

    # format blast result to have dict with all blast hit
    else:
        blast_result_dict = format_blast_output(out, gene_size_dict)
        return blast_result_dict

def format_blast_output(out, gene_size_dict):

    genome_result_dict = {}
    match = 1

    #checking length of single blast alignments against query. has to be 90% similar.
    if len(out.split('\n')) == 2:
        col = out.split('\t')

        query_length = query_length_gen(col)

        if query_length >= 0.8 * int(gene_size_dict[col[1]]):
            genome_result_dict[match] = out.strip('\n')
            match = match + 1
            return genome_result_dict

        else:
            genome_result_dict[match] = out.strip('\n') + '\t' + "Only blast hit represented < 90% of alignment in genome."
            match = match + 1
            return genome_result_dict


    #checking if allele blast across multiple contigs (unassembled region)
    if len(out.split('\n')) > 2:
        new_out_list = []

        total_length = 0
        sep_result = out.split('\n')[:-1]
        for res in sep_result:
            res_col = res.split('\t')

            query_length = query_length_gen(res_col)

            total_length = total_length + query_length

            new_out_list.append(res)

        new_out = '\t'.join(new_out_list)

        if total_length >= 0.8 * int(gene_size_dict[res_col[1]]):
            genome_result_dict[match] = new_out
            match = match + 1

            return genome_result_dict

        else:
            genome_result_dict[match] = new_out + '\t' + "Only blast hit represented < 90% of alignment in genome."
            match = match + 1
            return genome_result_dict

def query_length_gen(col):

    if int(col[9]) > int(col[8]):
        query_length = int(col[9]) - int(col[8]) + 1
    else:
        query_length = int(col[8]) - int(col[9]) + 1

    return query_length

def write_out(results_dict, args):
    print("Writing Out Results")

    #column names
    headers_list = ["Accession", "query", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

    #create output file and fill
    with open(args.output_folder + '/blast_results.csv','w') as out:

        #for each database type
        for key, value in results_dict.items():
            out.write(key + " Results" + '\n')

            #write out column headers
            for col_head in headers_list:
                out.write(col_head + ',')
            out.write('\n')

            #write out blast hits for each genome
            for genome in results_dict[key]:
                for match in results_dict[key][genome]:
                    col = results_dict[key][genome][match].split('\t')
                    out.write(genome + ',')
                    for cell in col:
                        out.write(cell + ',')
                    out.write('\n')
            out.write('\n')

def parseargs(set_wd):

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-db", "--database", nargs='+',
                        help="Choose which database to analyse genomes. ctx = 0, serogroup = 1, serotype = 2, tcp = 3, all databases = 4.\
                             For a combination use comman separated. Eg. 0,1,3.")
    parser.add_argument("-dir", "--strains_directory", required=True,
                        help="A directory of strains to analyse.")
    parser.add_argument("-r", "--blast_return", default=1, type=int,
                        help="1: return only top blast hits, 2: return all blast hits")
    parser.add_argument("-o", "--output_folder", default=set_wd + "/output/",
                        help="Output folder to save if not specified.)")


    args = parser.parse_args()

    return args

def main():

    set_wd = path.dirname(path.abspath(__file__))
    args = parseargs(set_wd)

    databases = choose_database(set_wd, args)
    results_dict = iterate_databases(databases, set_wd, args)
    write_out(results_dict, args)

if __name__ == '__main__':
    main()

#TODO fix serogroup output