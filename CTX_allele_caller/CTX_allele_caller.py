import argparse
import  glob
from os import path
import sys
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import multiprocessing
from time import sleep as sl

def run_single_genomes(set_wd, args):
    print("Processing genome :")

    #return single strain results as dict with number of blast hits
    results_dict = {}
    strain = args.input_genome
    accession = strain.split('/')[-1].split('.')[0]
    print(accession)
    results_dict[accession] = blast_input_against_db(strain, set_wd)

    #write out single strain
    write_out(results_dict, args)

def run_multiple_genomes(set_wd, args):
    print("Processing genomes")
    results_dict = {}
    for strain in glob.iglob(args.strains_directory + '*'):
        accession = strain.split('/')[-1].split('.')[0]
        print(accession)
        results_dict[accession] = blast_input_against_db(strain, set_wd, args)


    # write out single strain
    write_out(results_dict, args)

def blast_input_against_db(strain, set_wd, args):

    cpus = multiprocessing.cpu_count()
    blast_db = set_wd + "/blast_db/ctx_alleles_db"

    #only want top hit
    if args.blast_return == 1:
        blast_string = NcbiblastnCommandline(query=strain,db=blast_db, outfmt=6, perc_identity=50, num_threads=cpus, max_target_seqs=1)
        out, err = blast_string()

    #want all hits ##WARNING: many with low lengths
    elif args.blast_return == 2:
        blast_string = NcbiblastnCommandline(task="blastn", query=strain, db=blast_db, outfmt=6, perc_identity=50, num_threads=cpus)
        out, err = blast_string()

    #if not blast alignment found
    if out == '':
        blast_result_dict = {1:'No alignment found. Modify "-r 2" to see low quality hits.'}
        return blast_result_dict

    # format blast result to have dict with all blast hit
    else:
        blast_result_dict = format_blast_output(out)
        return blast_result_dict

def format_blast_output(out):

    ##format output into dict
    genome_result_dict = {}
    genome_result = out.split('\n')
    match = 1

    # for each returned blast hit
    for genome_pos in genome_result:
        if genome_pos != '':

            # write blast hit result to dictionary
            genome_result_dict[match] = genome_pos

            match = match + 1

    return genome_result_dict

def write_out(blast_results, args):
    print("Writing Out Results")

    #column names
    headers_list = ["query", "Accession", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

    #create output file and fill
    with open(args.output_folder + '/blast_results.csv','w') as out:
        #write out column headers
        for col_head in headers_list:
            out.write(col_head + ',')
        out.write('\n')

        #write out blast hits for each genome
        for key in blast_results:
            for match in blast_results[key]:
                col = blast_results[key][match].split('\t')
                out.write(key + ',')
                for cell in col:
                    out.write(cell + ',')
                out.write('\n')

def parseargs(set_wd):

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-i", "--input_genome",
                        help="Path to a single input genome (nucleotide) for analysis.")
    parser.add_argument("-dir", "--strains_directory",
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

    if args.input_genome:
        run_single_genomes(set_wd, args)

    elif args.strains_directory:
        run_multiple_genomes(set_wd, args)

    else:
        print("A single genome or list of genomes must be provided.")
        sys.exit()

if __name__ == '__main__':
    main()

#TODO upgrade snp comparison to just compare the SNP positions for each strain