#!/usr/bin/env python3
from time import sleep as sl

import sys
from os import environ
from Bio import SeqIO
import os
from os.path import commonprefix
import shutil
import argparse
# import psutil
import subprocess
from csv import reader
import multiprocessing
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from os.path import basename
from os import remove
import re
import itertools
import datetime
from copy import deepcopy
import time

from sys import stderr


import dis
from statistics import mean

"""
Dependencies:

Python 3
conda install psutil
conda install biopython
conda install Kraken 
(will need to DL minikraken DB)
export KRAKEN_DEFAULT_DB='/srv/scratch/lanlab/kraken_dir/minikraken_20141208'
conda install assembly-stats
conda install blast
conda install mafft
conda install -c bioconda sistr_cmd
conda install -c bioconda mlst

shovill (include 15cov version)

shovill dependencies
    skesa
    megahit (not used)
    velvet (not used)
    spades (not used)
    Lighter
    FLASH
    SAMtools >= 1.3
    BWA MEM
    MASH >= 2.0
    seqtk
    pigz
    Pilon (Java)
    Trimmomatic (Java)

conda:
Java
kraken



"""


# TODO return wrong input status (check that it is intact/valid fastq file) to the database (update db to say input incorrect)
# TODO also write back to db for contamination or genome quality


def run_mlst(ingenome):
    ##### Run 7 gene MLST program #####

    mlst_cmd = "mlst {}".format(ingenome)

    proc2 = subprocess.Popen(mlst_cmd, shell=True, stdout=subprocess.PIPE)

    mlst_result = proc2.stdout.read()
    mlst_result = mlst_result.decode('utf8')

    MGT1ST = mlst_result.split("\t")[2]

    return MGT1ST

######## ASSEMBLY TO ALLELES ########


def genome_to_alleles(query_genome, strain_name, args, mgt1st=""):
    """
    Takes assembly from assemblypipe and blasts against set of known alleles for all loci
    exact matches to existing alleles are called from perfect blast hits
    a locus that has no hits in the allele are called as 0

    :param query_genome: assembled genome path
    :param strain_name:
    :param args: args from main()
    :param mgt1st: 7 gene mlst included in allele file as a header

    :return: writes alleles,zero calls,7geneMLST to output file
    """
    start_time = time.time()

    # test_locus = "STM4037"
    stderr.write("Parsing inputs\n")
    script_path = sys.path[0]
    ref_alleles_in = args.refalleles  # known, intact, allele fasta file


    hsp_ident_thresh = float(args.hspident)  # scriptvariable blast identity to at least one other allele for each locus
    missing_limit = args.locusnlimit  # scriptvariable minimum allowable fraction of locus not lost (i.e. max 20% can be "N")
    wordsize = 18  # scriptvariable blast word size

    ####

    pref = args.outpath

    outdir = pref + strain_name

    outfile = outdir + "/" + strain_name + "_alleles.fasta"

    tempdir = outdir + "/tmp"

    if os.path.exists(outdir):
        pass
    else:
        os.mkdir(outdir)

    if os.path.exists(tempdir):
        shutil.rmtree(tempdir)
        os.mkdir(tempdir)
    else:
        os.mkdir(tempdir)

    qgenomeseq = SeqIO.parse(query_genome, "fasta")

    qgenome = {x.id: str(x.seq) for x in qgenomeseq}  # query genome as dictionary of {fastq_header:seq}

    locus_allowed_size,seqs,max_allele = get_allowed_locus_sizes(ref_alleles_in)

    locus_list = [x for x in locus_allowed_size.keys()]

    stderr.write("Running BLAST\n")

    # gets ref allele hits and blast hits against reference
    alleles_called_ref, ref_blast_hits, no_hits = ref_exact_blast(query_genome, ref_alleles_in, locus_allowed_size,
                                                                  tempdir)

    # print(no_hits)

    # print("Exact matches found: {}\n".format(len(alleles_called_ref.keys())))
    stderr.write("Processing partial BLAST hits\n")

    partial_loci = list(locus_list)  # list of all loci

    for i in alleles_called_ref:
        partial_loci.remove(i)  # if locus has exact match to existing remove from list

    # get hsps that match to loci but are not intact for partial_loci list
    # also dict of loci and reasons for uncallable loci in query genome
    partial_hsps, uncallable, tophitlocus = get_partial_match_query_region(ref_blast_hits, partial_loci, qgenome,
                                                                           hsp_ident_thresh,seqs)

    stderr.write("Reconstructing fragmented loci\n")

    # Try to rebuild each locus that has partial hsps matching it
    # returns reconstructed loci (with Ns) where possible
    reconstructed, uncallable = generate_query_allele_seqs(partial_hsps, query_genome, missing_limit,
                                                           wordsize, tophitlocus, qgenome, hsp_ident_thresh, uncallable)

    stderr.write("Writing outputs\n")
    write_outalleles(outfile, reconstructed, alleles_called_ref, uncallable, locus_list, mgt1st, no_hits,max_allele)

    elapsed_time = time.time() - start_time

    stderr.write("Allele calling completed in: {}\n".format(elapsed_time))

    now = datetime.datetime.now()
    timestamp = now.strftime("%H:%M:%S")

    # print("[" + timestamp + "] MGT fastq to alleles pipeline complete for strain: " + strain_name)


def ref_exact_blast(query_genome, ref_fasta, allele_sizes, tempdir):
    """
    runs blast and extracts exact hits to existing alleles
    data structure of parsed blast results:

    list of results (if multifasta input, one result per fasta seq)

    result attributes: >>>'alignments'<<<, 'application', 'blast_cutoff', 'database', 'database_length', 'database_letters', 'database_name', 'database_sequences', 'date', 'descriptions', 'dropoff_1st_pass', 'effective_database_length', 'effective_hsp_length', 'effective_query_length', 'effective_search_space', 'effective_search_space_used', 'expect', 'filter', 'frameshift', 'gap_penalties', 'gap_trigger', 'gap_x_dropoff', 'gap_x_dropoff_final', 'gapped', 'hsps_gapped', 'hsps_no_gap', 'hsps_prelim_gapped', 'hsps_prelim_gapped_attemped', 'ka_params', 'ka_params_gap', 'matrix', 'multiple_alignment', 'num_good_extends', 'num_hits', 'num_letters_in_database', 'num_seqs_better_e', 'num_sequences', 'num_sequences_in_database', 'posted_date', 'query', 'query_id', 'query_length', 'query_letters', 'reference', 'sc_match', 'sc_mismatch', 'threshold', 'version', 'window_size']
        alignment attributes: 'accession', 'hit_def', 'hit_id', >>>'hsps'<<<, 'length', 'title']
            hsp attributes: 'align_length', 'bits', 'expect', 'frame', 'gaps', 'identities', 'match', 'num_alignments', 'positives', 'query', 'query_end', 'query_start', 'sbjct', 'sbjct_end', 'sbjct_start', 'score', 'strand'

    :param query_genome: query genome fasta path
    :param ref_fasta: path to existing intact alleles file
    :param allele_sizes: dictionary of required allele sizes for each locus
    :param tempdir: directory path that is created to hold blast outputs/tmp files

    :return: dict of loci with exact hits, Bio.Blast.NCBIXML parsed blast results, list of loci with no blast hits in query genome
    """
    no_hits = list(allele_sizes.keys())
    # scriptvariable 15 blast word size
    # scriptvariable 1000000 culling limit
    # scriptvariable 90 blast identity limit

    bid = int(float(args.hspident)*100)

    blast_hits = run_blast(query_genome, ref_fasta, 15, 1000000, bid, tempdir)
    exact_list = []
    exact_dict = {}
    for result in blast_hits:
        contig_hit = result.query.split(" ")[0]
        for alignment in result.alignments:
            allele_hit = alignment.hit_def.split(":")[0]  # get locus name
            if allele_hit in allele_sizes:
                if allele_hit in no_hits:
                    no_hits.remove(allele_hit)  # remove locus from list of no hits
                allele_no = alignment.hit_def.split(":")[1]  # retreive allele number from full allele name
                for hsp in alignment.hsps:
                    if int(hsp.identities) == int(hsp.align_length) and int(hsp.gaps) == 0 and int(hsp.align_length) in \
                            allele_sizes[allele_hit]:
                        #  if no gaps, matching bases = length of alignment & length is same as input dictionary add to exact
                        exact_list.append(allele_hit)
                        if allele_hit not in exact_dict:
                            exact_dict[allele_hit] = allele_no  # store exact hit in dict

    return exact_dict, blast_hits, no_hits


def run_blast(query_seq, locus_db, wordsize, culling, pident, tempdir):
    """

    :param query_seq: query sequence - can be multiple fasta seqs
    :param locus_db: blastdb path
    :param wordsize: blast word size
    :return: returns list of blast results
    """
    # scriptvariable max hsps 5
    # scriptvariable evalue 0.1

    cpus = multiprocessing.cpu_count()
    gene = basename(locus_db).replace(".fasta", "")
    tmp_out = tempdir + "/tmp_blast.xml"
    cline = NcbiblastnCommandline(
        query=query_seq,
        subject=locus_db,
        evalue=0.1,
        perc_identity=pident,
        out=tmp_out,
        outfmt=5,
        max_target_seqs=100000,
        max_hsps=5,
        word_size=wordsize,
        num_threads=cpus,
        task="blastn",
        culling_limit=culling)
    # best_hit_score_edge = 0.01,
    try:
        stdout, stderr = cline()
    except Exception as e:
        print(e)
        sys.exit()

    r_handle = open(tmp_out)

    blast_records = list(NCBIXML.parse(r_handle))

    remove(tmp_out)
    """
    blast_records structure: 
    list of results (if multifasta input, one result per fasta seq) 

    result attributes: >>>'alignments'<<<, 'application', 'blast_cutoff', 'database', 'database_length', 'database_letters', 'database_name', 'database_sequences', 'date', 'descriptions', 'dropoff_1st_pass', 'effective_database_length', 'effective_hsp_length', 'effective_query_length', 'effective_search_space', 'effective_search_space_used', 'expect', 'filter', 'frameshift', 'gap_penalties', 'gap_trigger', 'gap_x_dropoff', 'gap_x_dropoff_final', 'gapped', 'hsps_gapped', 'hsps_no_gap', 'hsps_prelim_gapped', 'hsps_prelim_gapped_attemped', 'ka_params', 'ka_params_gap', 'matrix', 'multiple_alignment', 'num_good_extends', 'num_hits', 'num_letters_in_database', 'num_seqs_better_e', 'num_sequences', 'num_sequences_in_database', 'posted_date', 'query', 'query_id', 'query_length', 'query_letters', 'reference', 'sc_match', 'sc_mismatch', 'threshold', 'version', 'window_size']
        alignment attributes: 'accession', 'hit_def', 'hit_id', >>>'hsps'<<<, 'length', 'title']
            hsp attributes: 'align_length', 'bits', 'expect', 'frame', 'gaps', 'identities', 'match', 'num_alignments', 'positives', 'query', 'query_end', 'query_start', 'sbjct', 'sbjct_end', 'sbjct_start', 'score', 'strand'
    """

    return blast_records


def get_partial_match_query_region(blast_results, partial_matches, qgenome, hsp_ident,seqs):
    """

    :param blast_results:
    :param partial_matches: loci without exact matches
    :param sizes: not used
    :param ref_alleles:
    :param qgenome:
    :param hsp_ident:
    :return:
    """

    '''

    :param blast_results: blast results
    :param partial_matches: loci without exact matches
    :return: partials dictionary {locus:list of tuples} tuple -> (hsp matching reference allele,query contig where match was found)
    Also writes fasta file of these hits for each locus to be used to blast all alleles of the locus
    '''
    test = "-"
    partials = {}

    no_call_reason = {}

    # Combine below loop with finding best hit for each loci - store current best combined bitscore alignment in dict
    # and also store those hsps in partials - if higher score is found then replace in partials.

    top_allele_hits = {}

    tophitlocus = {}

    for result in blast_results:
        for alignment in result.alignments:
            alignment_locus = alignment.hit_def.rsplit(":")[0]
            full_subjct = seqs[alignment_locus][alignment.hit_def]
            if alignment_locus ==test:
                for hsp in alignment.hsps:
                    # print("\n")
                    print(alignment.hit_def)
                    print(hsp.query)
                    print(hsp.align_length)
                    print(hsp.sbjct)
                    print(hsp.score)
                    print("\n")
            if alignment_locus in partial_matches:  # and alignment.hit_def == alignment_locus + ":1":

                # if hsp locus is in partial matches list and allele is "ref" then store this hsp as partial for locus

                if alignment_locus not in top_allele_hits:
                    top_allele_hits[alignment_locus] = 0.0

                currentAlignBitscore = 0
                for hsp in alignment.hsps:
                    if hsp_filter_ok(hsp, 100, hsp_ident):
                        currentAlignBitscore += float(hsp.score)
                if alignment_locus ==test:
                    print(alignment.hit_def, currentAlignBitscore)

                if currentAlignBitscore > top_allele_hits[alignment_locus]:
                    top_allele_hits[alignment_locus] = currentAlignBitscore
                    partials[alignment_locus] = []

                    for hsp in alignment.hsps:
                        if hsp_filter_ok(hsp, 100, hsp_ident):

                            tophitlocus[alignment_locus] = full_subjct
                            partials[alignment_locus].append((hsp, result.query, alignment.hit_def))




    rmlis = []

    partials2 = {}

    # deal with partial overlaps of 2 hsps
    for locus in partial_matches:
        hsplis = []
        c = 1
        ## if 2 or more hsps pass filter but overlap by more that 60% call as 0

        if locus not in tophitlocus:
            rmlis.append(locus)
            no_call_reason[locus] = "no_blast_hits"
            olcheck = "nohit"

        elif len(partials[locus]) > 1:

            olcheck = check_for_multiple_ol_partial_hsps(partials[locus], hsp_ident)

            if olcheck == "overlap":
                rmlis.append(locus)
                no_call_reason[locus] = "large_overlap"
            else:
                locuspartials = remove_hsps_entirely_within_others(partials[locus], 0, hsp_ident)
        else:
            locuspartials = partials[locus]
            olcheck = "no overlap"
        if olcheck == "no overlap":
            # TODONE sizechange dict of allele sizes for each locus then use top hit to give size instead of ref_alleles
            locuspartials2 = check_ends_for_snps(locuspartials, tophitlocus[locus], locus, qgenome)
            partials2[locus] = locuspartials

            if len(locuspartials2) == 0:
                rmlis.append(locus)
                no_call_reason[locus] = "possible_duplication"

    npartials = {}
    for locus in partials2:
        if locus not in rmlis:
            npartials[locus] = partials2[locus]

    return npartials, no_call_reason, tophitlocus


def mask_high_snp_regions(locus, recon_locus, ref_locus, window_size, snp_limit):
    """
    compare new allele seq to "ref" allele to identify and mask(with Ns) regions with high snp density
    which can be caused by BLAST errors where indels are present in query seq

    :param locus: Not currently used (useful for debug)
    :param recon_locus: sequence to be checked
    :param ref_locus: sequence of 'ref' locus
    :param window_size: size of rolling window to check SNP frequency within
    :param snp_limit: limit of number of SNPs within that window

    :return: input sequence to be checked with regions with elevated SNP counts masked if necessary
    """

    if len(recon_locus) != len(ref_locus):  # comparison with reference will break if lengths not the same
        return recon_locus

    halfwindow = int(window_size / 2)
    mutpos = []
    outlocus = list(str(recon_locus))  # Convert test sequence into list of letters
    for pos in range(len(recon_locus)):
        if recon_locus[pos] != ref_locus[pos] and recon_locus[pos] not in ["N", "-"] and ref_locus[pos] not in ["N",
                                                                                                                "-"]:
            mutpos.append(
                "X")  # if the same position doesn't match (i.e. a SNP) and that missmatch is not caused by an N or an indel
        else:
            mutpos.append("M")  # If ref and new seq are the same

    #  Get list of ranges allowing for window length over whole locus
    for x in range(halfwindow, len(ref_locus) - halfwindow):
        window = mutpos[x - halfwindow:x + halfwindow]  # get window list of Match(M) or SNP(X)
        if window.count("X") > int(snp_limit):
            # if num of SNPS greater than limit mask all positions in current window
            for pos in range(x - halfwindow, x + halfwindow):
                outlocus[pos] = "N"  # change to N
    outlocus = "".join(outlocus)

    # if outlocus.count("N") > (1 - float(0.8)) * len(ref_locus):
    #     print(locus)

    return outlocus


def generate_query_allele_seqs(partial_hsps, query_genome, missing_perc_cutoff, wordsize, tophitlocus,
                               qgenome, hsp_thresh, uncallable):
    """

    reconstruct allele sequences from hsps that partially cover the locus

    :param partial_hsps: each locus' partially matching hsps
    :param query_genome: path to query genome file
    :param alleles_sizes: allele sizes dict
    :param missing_perc_cutoff: max amount of locus allowed to be missing as fraction (i.e. 0.8)
    :param wordsize: blast word size
    :param ref_alleles: allele seqs derived from "reference" (all are allele 1)
    :param qgenome: query genome as dict
    :param hsp_thresh: BLAST identity threshold to pass hsp filter
    :param uncallable: dict - loci called 0 with reason as value

    :return: reconstructed alleles in dict and loci where reconstruction failed in uncallable dict
    """
    query_genome = SeqIO.parse(query_genome, "fasta")
    q_genome = {}
    calls = {}
    test = "-"
    full_allele = ""

    for s in query_genome:
        seq = "".join([x for x in s.seq if x in ['N','T','C','G','A']])
        q_genome[s.id] = str(seq)

    # check that number of identities in blast hits is at least X fraction of normal reference allele length
    # need to check that at least x of allele is covered
    hspcov = 0
    for locus in partial_hsps:
        hspls = partial_hsps[locus]
        reflen = len(tophitlocus[locus])

        frac_covered = get_combined_hsp_coverage_of_ref_allele(reflen, hspls, hsp_thresh)  # annotation in function

        if frac_covered < float(missing_perc_cutoff):
            uncallable[locus] = "unscorable_too_much_missing"
            hspcov += 1
        hspls = list(remove_hsps_entirely_within_others(hspls, locus, hsp_thresh))
        # removes small hsps within larger ones caused by partial hits to non-orthologous but related genes

        hspls = check_ends_for_snps(hspls, tophitlocus[locus], "", qgenome)  # annotation in function

        hspls2 = []
        if len(hspls) == 0:
            uncallable[locus] = "failed_filter"
        elif len(hspls) > 1:
            contigs = {}
            for tup in hspls:
                hsp = tup[0]
                # scriptvariable hsp min size
                if hsp_filter_ok(hsp, 30, hsp_thresh):
                    hspls2.append(tup)
                    contig = tup[1].split(" ")[0]
                    if contig not in contigs:
                        contigs[contig] = [tup]
                    else:
                        contigs[contig].append(tup)

            # if hits come from > 1 contig
            if len(contigs.keys()) > 1:

                message, full_allele = check_split_over_contigs(hspls2, query_genome, reflen,
                                                                locus)  # annotation in function
                if message == "inconsistent_overlap":
                    uncallable[locus] = message
                else:
                    calls[locus] = full_allele

            else:

                for c in contigs:  # there will only be one contig
                    fixed_mid = check_mid(c, contigs[c], q_genome, wordsize, reflen, locus,tophitlocus[locus])  # annotation in function

                    if fixed_mid == "possible_insertion":
                        uncallable[locus] = "possible_insertion"
                    elif fixed_mid == "mixed_orientation":
                        uncallable[locus] = "mixed_orientation"
                    else:
                        calls[locus] = fixed_mid

        else:  # if only one hsp matches

            hsp = hspls[0][0]
            queryid = hspls[0][1].split(" ")[0]

            seq = remove_indels_from_hsp(hsp).query

            calls[locus] = str(seq)


    missingperc = 0
    calls2 = dict(calls)
    for locus in calls:
        reflen = float(len(tophitlocus[locus]))
        if float(len(calls[locus])) > 1.25 * reflen:  # this is not currentlyused but if indels used later will be important
            uncallable[locus] = "unscorable_too_long"
            # print(locus,len(calls[locus]),1.25 * reflen)
        elif float(len(calls[locus])) < 0.75 * reflen:  # this is not currentlyused but if indels used later will be important
            # print(locus, len(calls[locus]), 0.75 * reflen)
            uncallable[locus] = "unscorable_too_short"

        elif calls[locus].count("N") > (1 - float(missing_perc_cutoff)) * reflen:
            uncallable[locus] = "unscorable_too_much_missing"
            missingperc += 1

        else:
            calls2[locus] = calls[locus]

    return calls2, uncallable


######## PER LOCUS ########

def get_combined_hsp_coverage_of_ref_allele(reflen, hspls, hsp_thresh):
    """

    :param reflen: length of intact allele
    :param hspls: list of hsps matching
    :param hsp_thresh: minimum amount of locus below which 0 is called (i.e 0.8 = >80% must be present)
    :return: fraction of intact locus covered by partial hsps
    """
    range_lis = []
    for tup in hspls:
        hsp = tup[0]
        if hsp_filter_ok(hsp, 30, hsp_thresh):
            if hsp.sbjct_start > hsp.sbjct_end:
                range_lis.append((hsp.sbjct_end, hsp.sbjct_start))
            else:
                range_lis.append((hsp.sbjct_start, hsp.sbjct_end))
    mrange_lis = merge_intervals(
        range_lis)  # merge intervals by overlapping regions (i.e. 100-180 + 160-220 -> 100-220)
    tot_covered = 0
    for hit in mrange_lis:
        tot_covered += hit[1] - hit[0]
    fraction_covered = float(tot_covered) / float(reflen)
    return fraction_covered


def merge_intervals(intervals):
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []

    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]

            if int(higher[0]) - 1 <= int(lower[1]):
                upper_bound = max(int(lower[1]), int(higher[1]))
                merged[-1] = (int(lower[0]), upper_bound)  # replace by merged interval
            else:
                merged.append(higher)
    return merged


def hsp_filter_ok(hsp, length, fraction_snps_diff):
    """
    Check whether hsp is long enough (> min length) + is similar enough (fraction snp diff)
    ignores Ns in query or subject
    :param hsp:
    :param length: min length of match
    :param fraction_snps_diff: i.e. 0.9 = max allowed snp fraction is 10%
    :return: Boolean
    """
    ncount = 0
    for pos in range(hsp.align_length - 1):
        if hsp.query[pos] == "N" or hsp.sbjct[pos] == "N":
            ncount += 1
    if hsp.align_length > length:
        if (float(hsp.identities + hsp.gaps + ncount) / hsp.align_length) > fraction_snps_diff:
            return True
    return False


def check_ends_for_snps(hsplis, full_subj, locus, qgenome):
    """
    BLAST will not match ends of loci where a SNP is in the final position, therefore need to check with below code
    and insert correct nucleotide (from the query genome) if a SNP is present
    :param hsplis: list of high scoring pairs (hsps) from blast
    :param full_subj: allele subject
    :param locus: not needed except for error reporting
    :param qgenome: query genome

    :return: list of modified hsps with locus ending snps included
    """
    nhspls = []
    for hsp in hsplis:
        contig = hsp[1]
        allele = hsp[2]

        thsp = hsp[0]
        nhsp = hsp[0]

        contig2 = contig.split(" ")[0]
        contigseq = qgenome[contig2]
        subseq = full_subj

        trial = 0
        if thsp.sbjct_start > thsp.sbjct_end:
            if thsp.sbjct_end == 2:
                if thsp.query_end < len(contigseq):
                    q_snp = contigseq[thsp.query_end]
                else:
                    q_snp = "N"  # if only final nucleotide is deleted from allele add N TODO INDEL record this as a del when dels are included for this and next 3 instances
                s_snp = reverse_complement(subseq[0])
                nhsp.align_length = nhsp.align_length + 1
                nhsp.match = nhsp.match + " "
                nhsp.query = nhsp.query + q_snp
                nhsp.sbjct = nhsp.sbjct + s_snp
                nhsp.sbjct_end = nhsp.sbjct_end - 1
                nhsp.query_end = nhsp.query_end + 1
            if thsp.sbjct_start == len(subseq) - 1:
                # trial = 1
                if thsp.query_start - 2 >= 0:
                    q_snp = contigseq[thsp.query_start - 2]
                else:
                    q_snp = "N"
                s_snp = reverse_complement(subseq[-1])
                nhsp.align_length = nhsp.align_length + 1
                nhsp.match = " " + nhsp.match
                nhsp.query = q_snp + nhsp.query
                nhsp.sbjct = s_snp + nhsp.sbjct
                nhsp.sbjct_start = nhsp.sbjct_start + 1
                nhsp.query_start = nhsp.query_start - 1
        else:
            if thsp.sbjct_start == 2:
                if thsp.query_start - 2 >= 0:
                    q_snp = contigseq[thsp.query_start - 2]
                else:
                    q_snp = "N"
                s_snp = subseq[0]
                nhsp.align_length = nhsp.align_length + 1
                nhsp.match = " " + nhsp.match
                nhsp.query = q_snp + nhsp.query
                nhsp.sbjct = s_snp + nhsp.sbjct
                nhsp.sbjct_start = nhsp.sbjct_start - 1
                nhsp.query_start = nhsp.query_start - 1
            if thsp.sbjct_end == len(subseq) - 1:
                try:
                    if thsp.query_end < len(contigseq):
                        q_snp = contigseq[thsp.query_end]
                    else:
                        q_snp = "N"
                    s_snp = subseq[-1]
                    nhsp.align_length = nhsp.align_length + 1
                    nhsp.match = nhsp.match + " "
                    nhsp.query = nhsp.query + q_snp
                    nhsp.sbjct = nhsp.sbjct + s_snp
                    nhsp.sbjct_end = nhsp.sbjct_end + 1
                    nhsp.query_end = nhsp.query_end + 1
                except:
                    print(locus, "query > subject")
        nhspls.append((nhsp, contig, allele))

    '''
    blast hits structure: 
    list of results

    result attributes: >>>'alignments'<<<, 'application', 'blast_cutoff', 'database', 'database_length', 'database_letters', 'database_name', 'database_sequences', 'date', 'descriptions', 'dropoff_1st_pass', 'effective_database_length', 'effective_hsp_length', 'effective_query_length', 'effective_search_space', 'effective_search_space_used', 'expect', 'filter', 'frameshift', 'gap_penalties', 'gap_trigger', 'gap_x_dropoff', 'gap_x_dropoff_final', 'gapped', 'hsps_gapped', 'hsps_no_gap', 'hsps_prelim_gapped', 'hsps_prelim_gapped_attemped', 'ka_params', 'ka_params_gap', 'matrix', 'multiple_alignment', 'num_good_extends', 'num_hits', 'num_letters_in_database', 'num_seqs_better_e', 'num_sequences', 'num_sequences_in_database', 'posted_date', 'query', 'query_id', 'query_length', 'query_letters', 'reference', 'sc_match', 'sc_mismatch', 'threshold', 'version', 'window_size']
        alignment attributes: 'accession', 'hit_def', 'hit_id', >>>'hsps'<<<, 'length', 'title']
            hsp attributes: 'align_length', 'bits', 'expect', 'frame', 'gaps', 'identities', 'match', 'num_alignments', 'positives', 'query', 'query_end', 'query_start', 'sbjct', 'sbjct_end', 'sbjct_start', 'score', 'strand'

    '''
    return nhspls


def check_split_over_contigs(hspls, query_genome, reflen, locus):
    """
    Check if matches are at ends of contigs (or are that last non N part of the contigs)
    if so get any overlap that may have occured

    :param contigs_dict: dict of {contig_name:hsp(s) for that contig}
    :param query_genome:
    :param reflen: required length of locus
    :return: outcome of reconstruction and sequence (if successful, if not then empty string returned)
    """

    # generate list of hsps as tuples of locus position match (start, end)
    range_lis = []
    range_d = {}
    for tup in hspls:
        hsp = tup[0]

        if hsp.sbjct_start > hsp.sbjct_end:
            sst = hsp.sbjct_end
            send = hsp.sbjct_start
        else:
            sst = hsp.sbjct_start
            send = hsp.sbjct_end
        t = (sst, send)
        range_lis.append(t)
        range_d[t] = tup

    # sort list to be in order of match start position
    sorted_range_list = sorted(range_lis, key=lambda tup: tup[0])

    # get first hsp
    hsp1 = range_d[sorted_range_list[0]][0]

    newseq = hsp1.query

    if hsp1.sbjct_start > hsp1.sbjct_end:
        newseq = reverse_complement(hsp1.query)

    for rang in range(len(sorted_range_list) - 1):
        cur_range_end = sorted_range_list[rang][1]  # position of current hsp end
        next_range_start = sorted_range_list[rang + 1][0]  # position of next hsp start
        cur_hsp = range_d[sorted_range_list[rang]][0]  # current hsp range
        next_hsp = range_d[sorted_range_list[rang + 1]][0]  # next hsp range
        # check for overlap in current and next hsp in list
        if cur_range_end > next_range_start - 1:
            # if there is an overlap of two hsps check that region of overlap is identical, if not call 0
            # if consistent ignore the first part of the next hsp that makes up the overlap and add to new seq
            overlap = cur_range_end - next_range_start

            cur_query = cur_hsp.query

            if cur_hsp.sbjct_start > cur_hsp.sbjct_end:
                cur_query = reverse_complement(cur_query)

            next_query = next_hsp.query

            if next_hsp.sbjct_start > next_hsp.sbjct_end:
                next_query = reverse_complement(next_query)

            end_cur = cur_query[-1 * (overlap + 1):]

            start_next = next_query[:overlap + 1]

            if end_cur == start_next:
                newseq += next_query[overlap + 1:]
            else:
                return "inconsistent_overlap", ""

            # check overlap for identity

        elif cur_range_end == next_range_start - 1:
            # if the hsps end and start on consecutive positions then just add next hsp to new seq
            add = next_hsp.query
            if next_hsp.sbjct_start > next_hsp.sbjct_end:
                add = reverse_complement(next_hsp.query)
            newseq += add
        elif cur_range_end < next_range_start - 1:
            # if there is a gap add Ns the size of the gap to newseq and then the next hsp
            # add Ns in gap
            missing_no = next_range_start - cur_range_end - 1
            nstring = "N" * missing_no
            newseq += nstring
            add_seq = next_hsp.query
            if next_hsp.sbjct_start > next_hsp.sbjct_end:
                add_seq = reverse_complement(next_hsp.query)
            newseq += add_seq
        else:
            print("\n\nPROBLEM SPLIT CONTIGS\n\n")
    if len(newseq) == reflen:
        return "reconstructed split", newseq
    else:
        # this case caused by BLAST not calling missmatches at end of match
        # use check_ends_split_contigs to correct
        st_hsp = range_d[sorted_range_list[0]]
        en_hsp = range_d[sorted_range_list[-1]]

        full_allele, added_start, added_end = check_ends_split_contigs(st_hsp, en_hsp, reflen, newseq, locus)

        return "reconstructed split w ends", full_allele


def check_mid(contig, hspls, q_genome, wordsize, reflen, locus,refseq):
    """
    where two or more partial hsps match same contig
    sequentially add them together with Ns filling gaps and removing identical ovelaps to create full length allele
    use word length ratio as cutoff for allowing non-N letters in restored gap sequence(i.e. if word length is 12. non-N sequence chunk between 2 hsps can be a max of 18 before possible insertion is called)
    :param contig:
    :param hspls:
    :param q_genome:
    :param wordsize:
    :param reflen:
    :param locus:
    :return:
    """
    test = "-"
    # check that all hsps that match a locus are oriented the same way on the query genome contig
    orient = check_all_orient(hspls)
    genomekey = ""
    size = 0
    for con in q_genome:
        if len(q_genome[con]) > size:
            genomekey = con
            size = len(q_genome[con])

    if orient == "mixed":
        # TODO work out what to do with mixed orientation, currently call as 0
        return "mixed_orientation"
    else:

        order = {}
        for tup in hspls:
            hsp = tup[0]
            order[hsp.query_start] = hsp

        sorted_hsps = sorted(map(int, order.keys()))

        full_allele = remove_indels_from_hsp(order[sorted_hsps[0]]).query

        qstart = order[sorted_hsps[0]].query_start
        qend = order[sorted_hsps[-1]].query_end

        if orient == "positive":
            full_allele = q_genome[genomekey][qstart-1:qend]

        elif orient == "negative":
            full_allele = reverse_complement(q_genome[genomekey][qstart-1:qend])

        return full_allele


def check_all_orient(hsp_list):
    """
    Check orientation of hsps in list relative to subject ( existing alleles)
    and returns oritentation if the same or "mixed" if both orientations are found in same list
    :param hsp_list:
    :return:
    """
    orient = []
    for tup in hsp_list:
        hsp = tup[0]
        if hsp.sbjct_end < hsp.sbjct_start:
            orient.append("negative")
        else:
            orient.append("positive")
    if len(list(set(orient))) > 1:
        return "mixed"
    else:
        return orient[0]


def check_matching_overlap(hsp1, hsp2, orient,t):

    # if t:
    #     print(hsp1.sbjct_end,hsp2.sbjct_start)

    if orient == "positive":
        ol_len = hsp1.query_end - hsp2.query_start
        if hsp1.query[-1 * ol_len:] != hsp2.query[:ol_len]:
            return "nomatch", ol_len
        else:
            return "match", ol_len
    elif orient == "negative":
        ol_len = hsp1.query_end - hsp2.query_start
        if hsp1.query[:ol_len] != hsp2.query[-1 * ol_len:]:
            return "nomatch", ol_len
        else:
            return "match", ol_len


def remove_indels_from_hsp(hsp):
    """

    :param hsp:
    :return:
    """
    hsp1 = deepcopy(hsp)
    orient = "+"
    unknown_allele = hsp1.query
    ref_allele = hsp1.sbjct
    matches = hsp1.match
    refstart = hsp1.sbjct_start
    refend = hsp1.sbjct_end

    # re-orient matches so positions are always in direction of reference allele
    if int(hsp1.sbjct_start) > int(hsp1.sbjct_end):
        orient = "-"
        unknown_allele = reverse_complement(unknown_allele)
        ref_allele = reverse_complement(ref_allele)
        matches = matches[::-1]
        refstart, refend = refend, refstart

    seq = unknown_allele.replace("-","")
    hsp1.query = seq
    return hsp1


def check_ends(contig, qstart, qend, sstart, send, reflen, alleleseq, locus):
    """
    Add Ns to start and/or end
    :param contig: :param alleleseq:
    :param qstart: not used with no del calls
    :param qend: not used with no del calls
    :param sstart: hsp locus position start
    :param send: hsp locus position end
    :param reflen: required length of locus
    :param alleleseq: sequence of allele
    :param locus: not used with no del calls
    :return: reconstructed allele with sequences added, seq added to start, seq added to end
    """

    new_qstart = 0
    new_qend = 0
    added_start = ''
    added_end = ''

    full_allele = ""

    # if hit in reverse orientation have to check opposite end of hit in query genome

    if int(sstart) > int(send):

        # new_qstart = int(qstart) - (reflen - int(sstart))
        # new_qend = int(qend) + int(send) - 1
        added_start = "N" * (int(reflen) - (int(sstart)))
        # TODO Indel: following commented out elifs check for identity of missing region; if not Ns then assumed to be truncated allele - however when ignoring dels we can just add Ns to all of these regions becuase we treat dels and N regions the same
        '''
        # added_st = contig[new_qstart:int(qstart) - 1]
        # if added_st.count("N") == len(added_st):
        #     added_start = reverse_complement(added_st)
        # else:
        #     added_start = "N"*(int(reflen) - int(sstart))
        #     print("len of added start", int(reflen) - int(sstart))
        '''
        added_end = "N" * (int(send) - 1)
        '''
        # added_e = contig[int(qend):new_qend - 1]
        # if added_e.count("N") == len(added_e):
        #     added_end = reverse_complement(added_e)
        # else:
        #     added_end = "N"*(int(send)-1)
        #     print("len of added end",int(send)-1)
        '''
        full_allele = added_end + alleleseq + added_start


    elif int(sstart) < int(send):
        added_start = "N" * (int(sstart) - 1)
        added_end = "N" * (reflen - (int(send)))
        full_allele = added_start + alleleseq + added_end
        '''
        # new_qstart = int(qstart) - int(sstart)
        # new_qend = int(qend) + (reflen - int(send) - 1)

        # added_st = contig[new_qstart:int(qstart) - 1]
        # if added_st.count("N") == len(added_st):
        #     added_start = added_st
        # else:
        #     added_start = "N"*int(sstart-1)

        # added_e = contig[int(qend):new_qend - 1]
        # if added_e.count("N") == len(added_e):
        #     added_end = added_e
        # else:
        #     added_end = "N"*(reflen - int(send))
        '''

    # if missing flanking region is not all Ns then use original query hit edge - allows for rearrangements/deletions where flanking regions are really gone

    return full_allele, added_start, added_end


def check_ends_split_contigs(hsp_start_tup, hsp_end_tup, reflen, alleleseq, locus):
    """
    functions almost the same as check_ends()
    :param hsp_start_tup:
    :param hsp_end_tup:
    :param reflen:
    :param alleleseq:
    :return:
    """
    """
    If in correct orientation then start of start and end of end hsp will work
    BUT if one or both are in different orientation then need to swap start and end
    """

    # if hsp_start_tup[0] < hsp_start_tup[1]:
    #     start_hsp = hsp_start_tup[0]
    # else:
    #     start_hsp = hsp_start_tup[1]
    #
    # if hsp_end_tup[0] < hsp_end_tup[1]:
    #     end_hsp = hsp_end_tup[0]
    # else:
    #     end_hsp = hsp_end_tup[1]

    start_contig = hsp_start_tup[1]
    end_contig = hsp_end_tup[1]
    start_hsp = hsp_start_tup[0]
    end_hsp = hsp_end_tup[0]

    # start hsp

    # new_qstart = 0
    # new_qend = 0
    # sstart = start_hsp.sbjct_start
    # send = start_hsp.sbjct_end

    if start_hsp.sbjct_start < start_hsp.sbjct_end:
        sstart = start_hsp.sbjct_start
        send = start_hsp.sbjct_end
    else:
        sstart = start_hsp.sbjct_end
        send = start_hsp.sbjct_start

    qstart = start_hsp.query_start
    qend = start_hsp.query_end
    added_start = ""

    if int(sstart) > int(send):
        added_start = "N" * (int(send))
        # TODO Indel: following commented out elifs check for identity of missing region; if not Ns then assumed to be truncated allele - however when ignoring dels we can just add Ns to all of these regions becuase we treat dels and N regions the same
        '''
        # new_qstart = int(qstart) - (reflen - int(sstart))
        # added_st = start_contig[new_qstart:int(qstart) - 1]
        # if added_st.count("N") == len(added_st):
        #     added_start = reverse_complement(added_st)
        # else:
        #     added_start = "N"*len(added_st)
        # print(locus, "p3", added_start)
        '''
    else:
        added_start = "N" * (int(sstart) - 1)
        '''
        # new_qstart = int(qstart) - int(sstart)
        # added_st = start_contig[new_qstart:int(qstart) - 1]
        # if added_st.count("N") == len(added_st):
        #     added_start = added_st
        # else:
        #     added_start = "N"*len(added_st)
        # print(locus, "p4", added_start)
        '''

    # end_hsp

    if end_hsp.sbjct_start < end_hsp.sbjct_end:
        sstart2 = end_hsp.sbjct_start
        send2 = end_hsp.sbjct_end
    else:
        sstart2 = end_hsp.sbjct_end
        send2 = end_hsp.sbjct_start

    # sstart2 = end_hsp.sbjct_start
    # send2 = end_hsp.sbjct_end
    # qstart2 = end_hsp.query_start
    # qend2 = end_hsp.query_end

    added_end = ""

    if int(sstart2) < int(send2):
        added_end = "N" * (reflen - (int(send2)))
        '''
        # new_qend = int(qend2) + int(send2) - 1
        # added_e = end_contig[int(qend2):new_qend - 1]
        # if added_e.count("N") == len(added_e):
        #     added_end = reverse_complement(added_e)
        # else:
        #     added_end = "N"*len(added_e)
        '''
    elif int(sstart2) > int(send2):
        added_end = "N" * (reflen - (int(sstart2)))
        '''
        # new_qend = int(qend2) + (reflen - int(send2) - 1)
        # added_e = end_contig[int(qend2):new_qend - 1]
        # if added_e.count("N") == len(added_e):
        #     added_end = added_e
        # else:
        #     added_end = "N"*len(added_e)
        '''

    full_allele = added_start + alleleseq + added_end

    return full_allele, added_start, added_end


def get_allowed_locus_sizes(locusfile):
    sizes = {}
    seqs = {}
    file = SeqIO.parse(locusfile, "fasta")
    max_allele = 0
    for allele in file:
        locus = allele.id.split(":")[0]
        all = int(allele.id.split(":")[1])
        if all > max_allele:
            max_allele = all
        allelelen = len(allele.seq)
        if locus not in sizes:
            sizes[locus] = [allelelen]
            seqs[locus] = {allele.id:allele.seq}
        else:
            seqs[locus][allele.id] = allele.seq
            if allelelen not in sizes[locus]:
                sizes[locus].append(allelelen)
    return sizes,seqs,max_allele


######## PARTIAL BLAST HIT PROCESSING ########


def remove_hsps_entirely_within_others(hspls, locus, hsp_thresh):
    """

    :param hspls:
    :param reflen:
    :return:
    """
    range_d = {}
    range_d_new = {}

    for tup in hspls:

        hsp = tup[0]
        if hsp_filter_ok(hsp, 30, hsp_thresh):
            sst = 0
            send = 0
            if hsp.sbjct_start > hsp.sbjct_end:
                sst = hsp.sbjct_end
                send = hsp.sbjct_start
            else:
                sst = hsp.sbjct_start
                send = hsp.sbjct_end
            range_d[(sst, send)] = tup

    # if new range is entirely within existing one ignore
    # if new range entirely encompases existing one, remove existing and add new
    # otherwise add range
    remove = []
    for existing_range in range_d:
        for test_range in range_d:
            exst = existing_range[0]
            exend = existing_range[1]
            testst = test_range[0]
            testend = test_range[1]
            if exst <= testst and exend >= testend and existing_range != test_range:
                remove.append(test_range)
    nhspls = []
    for rang in range_d:
        if rang not in remove:
            range_d_new[rang] = range_d[rang]
            nhspls.append(range_d[rang])

    if len(nhspls) > 1:
        # TODONE write section to remove middle hsps like this:[(1, 220), (134, 285), (169, 498)] - i.e. remove (134, 285)
        r2_rangelist = []
        for rang in range_d_new:
            r2_rangelist.append(rang)

        r2_rangelist = sorted(r2_rangelist, key=lambda tup: tup[0])

        overlaps = []
        for r1 in r2_rangelist:
            list1_remove_r1 = [x for x in r2_rangelist if x != r1]
            for r2 in list1_remove_r1:
                if r1[1] > r2[0] and r1[1] < r2[1] and r1[0] < r2[0]:
                    overlaps.append((r1, r2))
        remove = []
        for ol in overlaps:
            for rang in r2_rangelist:
                if rang not in ol:
                    if ol[0][0] < rang[0] and rang[1] < ol[1][1]:
                        remove.append(rang)
        outls = []
        testls = []
        for x in r2_rangelist:
            if x not in remove:
                outls.append(range_d[x])
                testls.append(x)

        return outls
    else:
        return nhspls


def check_for_multiple_ol_partial_hsps(hspls, hsp_thresh):
    """
    check for large overlaps of hsps matching one locus currently cutoff is 60% length of larger hsp
    :param hspls: list of hsps
    :param hsp_thresh: hsp snp limit/ identity minimum (i.e. 0.98)
    :return:
    """
    # scriptvariable - hsp overlap % - cutoff for overlap elimination of partial hsps is currently 0.6 - could change / make variable to adjust

    # get list of hsps with start end and
    ## if 2 or more hsps pass filter but overlap by more that 60% call as 0
    hspls = [x[0] for x in hspls]
    filterpass = []
    for hsp in hspls:

        if hsp_filter_ok(hsp, 100, hsp_thresh):
            filterpass.append(hsp)

    for a, b in itertools.combinations(filterpass, 2):
        if a.align_length >= b.align_length:
            longer = a
            shorter = b
        else:
            longer = b
            shorter = a
        longerst = longer.sbjct_start
        longeren = longer.sbjct_end
        if longer.sbjct_start > longer.sbjct_end:
            longerst = longer.sbjct_end
            longeren = longer.sbjct_start
        shorterst = shorter.sbjct_start
        shorteren = shorter.sbjct_end
        if shorter.sbjct_start > shorter.sbjct_end:
            shorterst = shorter.sbjct_end
            shorteren = shorter.sbjct_start

        if longerst <= shorterst <= longeren <= shorteren:
            if float(longeren - shorterst) >= float(0.6 * longer.align_length):
                return "overlap"
        elif shorterst <= longerst <= shorteren <= longeren:
            if float(shorteren - longerst) >= float(0.6 * longer.align_length):
                return "overlap"
        elif longerst <= shorterst <= shorteren <= longeren:
            if float(shorter.align_length) >= float(0.6 * longer.align_length):
                return "overlap"
    return "no overlap"


######## UTILS ########


def reverse_complement(dna):
    """

    :param dna:
    :return:
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', "-": "-", "N": "N"}

    dna = ''.join([x for x in dna if x in ['N','A','T','G','C','-']])

    out = ''.join([complement[base] for base in dna[::-1]])

    return out


def largest_nonn_strings(string):
    """

    :param string:
    :return:
    """
    matches = re.findall(r"[ATGC]*", string)
    mx = max([len(x) for x in matches])
    return mx


def check_zp(fq1):
    if ".gz" in fq1:
        return True
    return False


######## IO ########


def get_sizes_dict(inf):
    infile = open(inf, "r").read().splitlines()
    infile = [x.split("\t") for x in infile]
    outd = {x[0]: str((int(x[2]) - int(x[1])) + 1) for x in infile}
    return outd


def write_outalleles(outpath, reconstructed, ref, uncall, locuslist, mgt1st, no_hits,max_allele):
    # outf = open(outpath, "w")
    # outf.write(">{}:{}\n\n".format("7_gene_ST", mgt1st))
    call = 0
    new = {}
    missing = {}
    absent = 0
    if args.newno:
        if args.newno.isdigit():
            nallele = str(args.newno)
        else:
            sys.exit("new allele number must be integer")
    else:
        nallele = max_allele+1
    for locus in locuslist:
        if locus in ref:
            # outf.write(">{}:{}\n\n".format(locus, ref[locus]))
            call += 1
        elif locus in uncall:
            # outf.write(">{}:0_{}\n\n".format(locus, uncall[locus]))
            if uncall[locus] not in missing:
                missing[uncall[locus]] = 1
            else:
                missing[uncall[locus]] += 1
        elif locus in reconstructed:
            print(">{}:{}\n{}".format(locus,nallele, reconstructed[locus]))
            # outf.write(">{}:{}\n{}\n".format(locus,nallele, reconstructed[locus]))
            if "N" in reconstructed[locus]:
                if "partial" not in new:
                    new["partial"] = 1
                else:
                    new["partial"] += 1
            else:
                if "intact" not in new:
                    new["intact"] = 1
                else:
                    new["intact"] += 1

        elif locus in no_hits:
            # outf.write(">{}:0_{}\n\n".format(locus, "no_blast_hits"))
            absent += 1
        else:
            # outf.write(">{}:0_{}\n\n".format(locus, "no_result"))
            # print("MISSING: ", locus)
            absent += 1
    stderr.write("called:{}\n".format(call))
    for i in new:
        stderr.write("{}: {}\n".format(i, new[i]))
    for i in missing:
        stderr.write("{}: {}\n".format(i, missing[i]))
    stderr.write("absent:{}\n\n".format(absent))
    # outf.close()


######## ARGUMENTS/HELP ########


if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("ingenome", help="Input genome")
    parser.add_argument("refalleles", help="File path to MGT reference alleles file")
    parser.add_argument("outpath", help="Path to ouput file name")

    parser.add_argument("-s", "--species", help="String to find in kraken species confirmation test",
                        default="Salmonella enterica")
    parser.add_argument("--no_serotyping", help="Do not run Serotyping of Salmonella using SISTR (ON by default)")
    parser.add_argument("-y", "--serotype", help="Serotype to match in SISTR, semicolon separated",
                        default="Typhimurium;I 4,[5],12:i:-")
    parser.add_argument("-t", "--threads", help="number of computing threads",
                        default="4")
    parser.add_argument("-m", "--memory", help="memory available in GB",
                        default="8")
    parser.add_argument("-n", "--newno", help="allele number to give to new alleles",type=int)
    parser.add_argument("-f", "--force", help="overwrite output files with same strain name?", action='store_true')
    parser.add_argument("--min_largest_contig",
                        help="Assembly quality filter: minimum allowable length of the largest contig in the assembly in bp (default for salmonella)",
                        default=60000)
    parser.add_argument("--max_contig_no",
                        help="Assembly quality filter: maximum allowable number of contigs allowed for assembly (default for salmonella)",
                        default=700)
    parser.add_argument("--genome_min",
                        help="Assembly quality filter: minimum allowable total assembly length in bp (default for salmonella)",
                        default=4500000)
    parser.add_argument("--genome_max",
                        help="Assembly quality filter: maximum allowable total assembly length in bp (default for salmonella)",
                        default=5500000)
    parser.add_argument("--n50_min",
                        help="Assembly quality filter: minimum allowable n50 value in bp (default for salmonella)",
                        default=20000)
    parser.add_argument("--kraken_db",
                        help="path for kraken db (if KRAKEN_DEFAULT_DB variable has already been set then ignore)",
                        default="")
    parser.add_argument("--hspident",
                        help="BLAST percentage identity needed for hsp to be returned",
                        default=0.98)
    parser.add_argument("--locusnlimit",
                        help="minimum proportion of the locus length that must be present (not masked with Ns)",
                        default=0.80)
    parser.add_argument("--snpwindow",
                        help="Size of sliding window to screen for overly dense SNPs",
                        default=40)
    parser.add_argument("--densitylim",
                        help="maximum number of SNPs allowed to be present in window before window is masked",
                        default=4)
    parser.add_argument("--refsize",
                        help="Approx size of genome for shovill input in megabases i.e. 5.0 or 2.9", type=float,
                        default=5.0)


    args = parser.parse_args()
    name = os.path.basename(args.ingenome).strip(".fasta")

    # mlst7gene = run_mlst(args.ingenome)

    genome_to_alleles(args.ingenome, name, args)
