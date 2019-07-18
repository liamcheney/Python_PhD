from time import sleep as sl
import os
import shutil
import itertools
import sys
from os import path
from Bio import SeqIO
import re
from copy import deepcopy
import time

from allele_MLST_caller import IO,wrappers,utils


sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) ) )

def hsp_filter_ok(hsp,length,fraction_snps_diff):
    ncount = 0
    for pos in range(hsp.align_length-1):
        if hsp.query[pos] == "N" or hsp.sbjct[pos] == "N":
            ncount +=1
    if hsp.align_length > length:
        if (float(hsp.identities+hsp.gaps+ncount) / hsp.align_length) > fraction_snps_diff:
            return True
    return False

def make_intact_allele_fasta(scheme_fastas,dbname,all_alleles,allele_subset_dict,scheme):
    """
    :param scheme_fastas: fasta file paths for scheme of interest
    :param dbname: output fasta path
    :param all_alleles
    :param allele_subset_dict
    :param scheme
    :return:fasta file with all intact alleles derived from fasta files of loci making up a scheme
    also return dict of {locus:(list of intact alleles,list of all alleles)}
    """


    alleles_dict = {}
    ref_alleles = {}
    allele_seqs = {}
    outseqs = []

    c=0
    all = 0
    for fname in scheme_fastas:
        t=0

        infile = SeqIO.parse(fname,"fasta")
        intact = []
        allalleles = []
        locus = fname.split("/")[-1].replace(".fasta","")

        for allele in infile:
            # if t ==1:print(allele,all_alleles)
            allele_seqs[allele.id] = str(allele.seq)

            all_alleles2 = bool(all_alleles)
            if not all_alleles:
                # print(allele_subset_dict[scheme][locus])
                allele_nos = list(set(allele_subset_dict[scheme][locus]))

                allele_nos = [x.split(":")[-1] for x in allele_nos]
                if allele_nos[0] == "0" and len(allele_nos) == 1:
                    all_alleles2 = True
            if "-" not in allele.id and allele.id != "0":
                #outseqs.append(allele)
                intactno = allele.id.split(":")[-1]
                if all_alleles2:
                    outseqs.append(allele)
                else:
                    if allele.id in allele_subset_dict[scheme][locus]:
                        outseqs.append(allele)

                intact.append(intactno)
                allalleles.append(intactno)
                c+=1
                all+=1
            else:
                alleleno = allele.id.split(":")[-1]
                # intac = allele.id.split(":")[-1].split("_")[0][1:]
                # intact.append(intac)
                allalleles.append(alleleno)
                all+=1
            if allele.id.split(":")[-1] == "1":
                ref_alleles[locus] = str(allele.seq)
                # outseqs.append(allele)
                # intact.append("1")
                # allalleles.append("1")
        alleles_dict[locus] = (intact,allalleles)

    #FIXME srites out temporary fasta file for BLASTING against for each scheme
    SeqIO.write(outseqs,dbname,"fasta")

    return alleles_dict, ref_alleles, allele_seqs


def call_indels(locus,hsp,indels_calls):
    unknown_allele = hsp.query
    ref_allele = hsp.sbjct
    matches = hsp.match
    refstart = hsp.sbjct_start
    refend = hsp.sbjct_end

    if int(hsp.sbjct_start) > int(hsp.sbjct_end):
        unknown_allele = utils.reverse_complement(unknown_allele)
        ref_allele = utils.reverse_complement(ref_allele)
        matches = matches[::-1]
        refstart, refend = refend, refstart

    intact_nucs = ["A", "T", "G", "C", "a", "t", "c", "g"]
    for i in range(len(unknown_allele)):
        if matches[i] == " ":
            if unknown_allele[i] == "-" and ref_allele[i] in intact_nucs:
                pos = i - ref_allele[:i].count("-")
                pos_rel_to_ref_allele = pos + refstart - 1
                d = (ref_allele[i], str(pos_rel_to_ref_allele), "del")
                if locus not in indels_calls:
                    indels_calls[locus] = []
                indels_calls[locus].append(d)
                # TODO INDEL REMOVE Below replaces deletion in unknown allele with reference nucleotides to IGNORE DELETIONS
            elif unknown_allele[i] in intact_nucs and ref_allele[i] == "-":
                pos = i - ref_allele[:i].count("-")
                pos_rel_to_ref_allele = pos + refstart - 1
                ins = ("ins", str(pos_rel_to_ref_allele), unknown_allele[i])
                if locus not in indels_calls:
                    indels_calls[locus] = []
                indels_calls[locus].append(ins)
                # TODO INDEL REMOVE No addition of unknown allele seq if insertion in new allele
    return indels_calls


def get_exactmatches(blasthits,partial_hit,allele_sizes,hit_locations,indels,ref_called_alleles,hsp_thresh):
    """

    :param blasthits: parsed blast output as list of qury seq results
    :param partial_hit:list of loci, when an exact hit is found it is removed from list - resulting in a list of loci still to be examined
    :param allele_sizes: allele lengths of reference alleles
    :param hit_locations
    :return: alleles that match exactly for length with no gaps and 100% identity - one for each locus - alleles with Ns are ignored
    """
    #TODO CKHECK IF  NECESSARY make so if perfect hit to 2 alleles (because of 2 alleles w diff sizes) store both and select hit with largest score - also check that 2 hits are in same locus otherwise 2 diff alleles in one strain - call as 0
    exacthits = 0
    totalhits = 0
    perfect_hit_ls = {}

    for result in blasthits:
        contig_hit = result.query.split(" ")[0]
        for alignment in result.alignments:
            allele_hit = alignment.hit_def
            locus = alignment.hit_def.split(":")[0]
            if locus not in ref_called_alleles:
                for hsp in alignment.hsps:
                    # if "STMMW_44761" in alignment.title:
                    #     if hsp.align_length > 50:
                    #         locus = alignment.hit_def.split(":")[0]
                    #         print(locus, allele_sizes[alignment.hit_def])
                    #         print("test1", (int(hsp.identities) + int(hsp.gaps)), int(hsp.align_length))
                    #         print("test2", int(hsp.gaps))
                    #         print("test3", (len(hsp.sbjct) - hsp.sbjct.count("-")), int(allele_sizes[alignment.hit_def]))
                    #         print_blasthsp(result, alignment, hsp, "")
                    #         sl(0.3)
                    if hsp_filter_ok(hsp, 30, hsp_thresh):

                        if int(hsp.identities) == int(hsp.align_length) and int(hsp.gaps) == 0 and int(hsp.align_length) == int(allele_sizes[alignment.hit_def]):
                            # if testlocus in alignment.title:
                            #     print("perfect_hit")
                            #     locus = alignment.hit_def.split("-")[0]
                            #     print(locus, allele_sizes[alignment.hit_def])
                            #     print_blasthsp(result, alignment, hsp, "")
                            perfect_hit_allele = alignment.title.split(" ")[-1]
                            perfect_hit_locus = perfect_hit_allele.split(":")[0]

                            if perfect_hit_locus not in perfect_hit_ls:
                                perfect_hit_ls[perfect_hit_locus] = [(hsp,allele_hit,contig_hit,perfect_hit_allele)]
                            else:
                                perfect_hit_ls[perfect_hit_locus].append((hsp,allele_hit,contig_hit,perfect_hit_allele))


                            # print(perfect_hit_locus)
                            # if perfect_hit_locus not in partial_hit:
                            #     if perfect_hit[perfect_hit_locus] != perfect_hit_allele:
                            #         perfect_hit[perfect_hit_locus] = "0"
                            #         #if second or more hits in genome match first perfect hit then leave as first hit call
                            #         #Multiple perfect hits in genome with different alleles - don't call ie.e call as allele 0
                            # else:
                            #     perfect_hit[perfect_hit_locus] = perfect_hit_allele
                            #     partial_hit.remove(perfect_hit_locus)
                            exacthits +=1
                            totalhits += 1
                            hit_locations[result.query.split(" ")[0]].append((hsp.query_start,hsp.query_end,1.00))
                        elif (int(hsp.identities)+int(hsp.gaps)) == int(hsp.align_length) and int(hsp.gaps) > 0 and (len(hsp.sbjct)-hsp.sbjct.count("-")) == int(allele_sizes[alignment.hit_def]):
                            # TODO INDEL REMOVE this elif ignores gaps that are not caused by Ns: effectively if allele matches based on SNPs then matched - ignores indels
                            #TODO record dels that are ignored here
                            # if testlocus in alignment.title:
                            #     print("passed check")
                            #     locus = alignment.hit_def.split("-")[0]
                            #     print(locus, allele_sizes[alignment.hit_def])
                            #     print_blasthsp(result, alignment, hsp, "")
                            if "N" not in hsp.query:
                                perfect_hit_allele = alignment.title.split(" ")[-1]
                                perfect_hit_locus = perfect_hit_allele.split(":")[0]


                                if perfect_hit_locus not in perfect_hit_ls:
                                    perfect_hit_ls[perfect_hit_locus] = [(hsp, allele_hit, contig_hit,perfect_hit_allele)]
                                else:
                                    perfect_hit_ls[perfect_hit_locus].append((hsp, allele_hit, contig_hit,perfect_hit_allele))
                                # print(perfect_hit_locus)
                                # if perfect_hit_locus not in partial_hit:
                                #     if perfect_hit[perfect_hit_locus] != perfect_hit_allele:
                                #         perfect_hit[perfect_hit_locus] = "0"
                                #         # if second or more hits in genome match first perfect hit then leave as first hit call
                                #         # Multiple perfect hits in genome with different alleles - don't call ie.e call as allele 0
                                # else:
                                #     perfect_hit[perfect_hit_locus] = perfect_hit_allele
                                #     partial_hit.remove(perfect_hit_locus)
                                exacthits += 1
                                totalhits += 1
                                hit_locations[result.query.split(" ")[0]].append((hsp.query_start, hsp.query_end, 1.00))
                            else:
                                totalhits += 1
                        elif hsp.identities > (hsp.align_length*0.5) and hsp.identities < hsp.align_length:
                            # if testlocus in alignment.title:
                            #     print("passed_to_partial")
                            #     locus = alignment.hit_def.split("-")[0]
                            #     print(locus, allele_sizes[alignment.hit_def])
                            #     print_blasthsp(result, alignment, hsp, "")
                            totalhits +=1
    perfect_hit = {}

    for locus in ref_called_alleles:
        if locus in partial_hit:
            perfect_hit[locus] = locus + ":"+ ref_called_alleles[locus]
            partial_hit.remove(locus)

    for locus in perfect_hit_ls:
        if locus in allele_sizes:
            if len(perfect_hit_ls[locus]) == 1:
                indels = call_indels(locus,perfect_hit_ls[locus][0][0],indels)
                perfect_hit[locus] = perfect_hit_ls[locus][0][3]
                partial_hit.remove(locus)
            elif len(perfect_hit_ls[locus]) > 1:
                ls_of_contigs = [x[2] for x in perfect_hit_ls[locus]]
                num_contigs = len(set(ls_of_contigs))
                if num_contigs > 1:
                    partial_hit.remove(locus) ## by removing from partial hit and not recording perfect hit this locus is assigned 0 by default
                    print("perfect_hits 2 contigs",locus)
                else:
                    for hit in perfect_hit_ls[locus]:
                        # check if hit shares allele start and end nos if none do then do nothing (will remain in partial hits and be assigned new allele)
                        hsp = hit[0]
                        st = hsp.sbjct_start
                        en = hsp.sbjct_end
                        if hsp.sbjct_start > hsp.sbjct_end:
                            st = hsp.sbjct_end
                            en = hsp.sbjct_start
                        if st == 1 and en == allele_sizes[locus]:
                            if locus in perfect_hit:
                                if allele_sizes[locus+":"+hit[3]] > allele_sizes[locus+":"+perfect_hit[locus]]:
                                    ## if two hits on same contig likely overlapping hits with different sizes - pick larger hit
                                    perfect_hit[locus]=hit[3]
                                    indels = call_indels(locus, hit[0], indels)
                            else:
                                perfect_hit[locus] = hit[3]
                                indels = call_indels(locus, hit[0], indels)
                                partial_hit.remove(locus)


    # print("exact hits: " + str(exacthits) + ", total hits: " + str(totalhits) + ", remaining loci to check: " + str(len(partial_hit)))
    # if testlocus in partial_hit:
    #     print("error is after exact match")
    # if "STMMW_40231" in partial_hit:
    #     print("partial")

    return partial_hit, perfect_hit,hit_locations,indels


def check_for_multiple_ol_partial_hsps(hspls,hsp_thresh):
    ##TODO NOTE cutoff for overlap elimination of partial hsps is currently 0.6 - could change / make variable to adjust
    # may be an issue with pairwise comparison callin g

    #get list of hsps with start end and
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
            if float(longeren - shorterst) >= float(0.6*longer.align_length):
                return "overlap"
        elif shorterst <= longerst <= shorteren <= longeren:
            if float(shorteren - longerst) >= float(0.6 * longer.align_length):
                return "overlap"
        elif longerst <= shorterst <= shorteren <= longeren:
            if float(shorter.align_length) >= float(0.6 * longer.align_length):
                return "overlap"
    return "no overlap"


def remove_hsps_entirely_within_others(hspls,locus,hsp_thresh):
    """

    :param hspls:
    :param reflen:
    :return:
    """
    range_d = {}
    range_d_new = {}

    for tup in hspls:

        hsp = tup[0]
        if hsp_filter_ok(hsp,30,hsp_thresh):
            # print(tup[1])
            # print('alignment length:', hsp.align_length)
            # print('identities:', hsp.identities)
            # print('gaps:', hsp.gaps)
            # print('ref_allele_len', str(reflen))
            # print('query start', hsp.query_start)
            # print('query end: ', hsp.query_end)
            # print('subject start', hsp.sbjct_start)
            # print('subject end: ', hsp.sbjct_end)
            # print(hsp.query)
            # print(hsp.match)
            # print(hsp.sbjct, "\n")
            sst = 0
            send = 0
            if hsp.sbjct_start > hsp.sbjct_end:
                sst = hsp.sbjct_end
                send = hsp.sbjct_start
            else:
                sst = hsp.sbjct_start
                send = hsp.sbjct_end
            range_d[(sst,send)] = tup

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
        #TODONE write section to remove middle hsps like this:[(1, 220), (134, 285), (169, 498)] - i.e. remove (134, 285)
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


def remove_indels_from_hsp(hsp):
    """

    :param hsp:
    :return:
    """
    hsp1 = hsp
    orient = "+"
    unknown_allele = hsp1.query
    ref_allele = hsp1.sbjct
    matches = hsp1.match
    refstart = hsp1.sbjct_start
    refend = hsp1.sbjct_end

    # re-orient matches so positions are always in direction of reference allele
    if int(hsp1.sbjct_start) > int(hsp1.sbjct_end):
        orient = "-"
        unknown_allele = utils.reverse_complement(unknown_allele)
        ref_allele = utils.reverse_complement(ref_allele)
        matches = matches[::-1]
        refstart, refend = refend, refstart

    seq = ""
    intact_nucs = ["A", "T", "G", "C", "a", "t", "c", "g"]
    for i in range(len(unknown_allele)):
        if matches[i] == " ":
            if unknown_allele[i] in intact_nucs and ref_allele[i] in intact_nucs:
                seq += unknown_allele[i]
            elif unknown_allele[i] == "-" and ref_allele[i] in intact_nucs:
                # TODO INDEL REMOVE Below replaces deletion in unknown allele with reference nucleotides to IGNORE DELETIONS
                seq += ref_allele[i]
            elif unknown_allele[i] in intact_nucs and ref_allele[i] == "-":
                continue
                # TODO INDEL REMOVE No addition of unknown allele seq if insertion in new allele
                # seq+=unknown_allele[i]
            elif unknown_allele[i] == "N":
                seq += "N"
        else:
            seq += unknown_allele[i]
    return seq


def check_ends_for_snps(hsplis,full_subj,query,locus,qgenome):
    #TODONE fix so that correct nucleotides are pulled across - also check that indels don't mess with this calc
    nhspls = []
    for hsp in hsplis:
        contig = hsp[1]
        allele = hsp[2]
        # print(contig,allele,query)
        thsp = hsp[0]
        nhsp = hsp[0]
        if query == "":
            contig2 = contig.split(" ")[0]
            contigseq = qgenome[contig2]
            subseq = full_subj
        else:

            # print(full_subj)
            contig = contig.replace("reconstructed_allele","")
            contigseq = full_subj[contig]
            subseq = query
        trial = 0
        if thsp.sbjct_start > thsp.sbjct_end:
            if thsp.sbjct_end == 2:
                if thsp.query_end < len(contigseq):
                    q_snp = contigseq[thsp.query_end] # TODONE had error here see 23-2-18 notes - 2-3-18 DONE
                else:
                    q_snp = "N" # if only final nucleotide is deleted from allele add N TODO indel record this as a del when dels are included for this and next 3 instances
                s_snp = utils.reverse_complement(subseq[0])
                nhsp.align_length = nhsp.align_length + 1
                nhsp.match = nhsp.match + " "
                nhsp.query = nhsp.query + q_snp
                nhsp.sbjct = nhsp.sbjct + s_snp
                nhsp.sbjct_end = nhsp.sbjct_end-1
                nhsp.query_end = nhsp.query_end+1
            if thsp.sbjct_start == len(subseq)-1:
                # trial = 1
                if thsp.query_start-2 >= 0:
                    q_snp = contigseq[thsp.query_start-2]
                else:
                    q_snp = "N"
                s_snp = utils.reverse_complement(subseq[-1])
                nhsp.align_length = nhsp.align_length + 1
                nhsp.match = " " + nhsp.match
                nhsp.query = q_snp + nhsp.query
                nhsp.sbjct = s_snp + nhsp.sbjct
                nhsp.sbjct_start = nhsp.sbjct_start+1
                nhsp.query_start = nhsp.query_start-1
        else:
            if thsp.sbjct_start == 2:
                if thsp.query_start-2 >= 0:
                    q_snp = contigseq[thsp.query_start-2]
                else:
                    q_snp = "N"
                s_snp = subseq[0]
                nhsp.align_length = nhsp.align_length + 1
                nhsp.match = " " + nhsp.match
                nhsp.query = q_snp + nhsp.query
                nhsp.sbjct = s_snp + nhsp.sbjct
                nhsp.sbjct_start = nhsp.sbjct_start-1
                nhsp.query_start = nhsp.query_start-1
            if thsp.sbjct_end == len(subseq)-1:
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
                    nhsp.sbjct_end = nhsp.sbjct_end+1
                    nhsp.query_end = nhsp.query_end+1
                except:
                    print(locus,"query > subject")
                    # print(locus)
                    # print_blasthsp("", "", nhsp, "")
                    # print(subseq)
                    # print(contigseq)
        nhspls.append((nhsp,contig,allele))



    '''
    blast hits structure: 
    list of results

    result attributes: >>>'alignments'<<<, 'application', 'blast_cutoff', 'database', 'database_length', 'database_letters', 'database_name', 'database_sequences', 'date', 'descriptions', 'dropoff_1st_pass', 'effective_database_length', 'effective_hsp_length', 'effective_query_length', 'effective_search_space', 'effective_search_space_used', 'expect', 'filter', 'frameshift', 'gap_penalties', 'gap_trigger', 'gap_x_dropoff', 'gap_x_dropoff_final', 'gapped', 'hsps_gapped', 'hsps_no_gap', 'hsps_prelim_gapped', 'hsps_prelim_gapped_attemped', 'ka_params', 'ka_params_gap', 'matrix', 'multiple_alignment', 'num_good_extends', 'num_hits', 'num_letters_in_database', 'num_seqs_better_e', 'num_sequences', 'num_sequences_in_database', 'posted_date', 'query', 'query_id', 'query_length', 'query_letters', 'reference', 'sc_match', 'sc_mismatch', 'threshold', 'version', 'window_size']
        alignment attributes: 'accession', 'hit_def', 'hit_id', >>>'hsps'<<<, 'length', 'title']
            hsp attributes: 'align_length', 'bits', 'expect', 'frame', 'gaps', 'identities', 'match', 'num_alignments', 'positives', 'query', 'query_end', 'query_start', 'sbjct', 'sbjct_end', 'sbjct_start', 'score', 'strand'

    '''
    return nhspls


def get_partial_match_query_region(blast_results,partial_matches,sizes,ref_alleles,qgenome,hsp_ident):
    """

    :param blast_results: 1st round blast results
    :param partial_matches: loci without exact matches
    :return: partials dictionary {locus:list of tuples} tuple -> (hsp matching reference allele,query contig where match was found)
    Also writes fasta file of these hits for each locus to be used to blast all alleles of the locus
    """

    # print(partial_matches)
    partials = {}
    for result in blast_results:
        for alignment in result.alignments:
            alignment_locus = alignment.hit_def.rsplit(":")[0]
            if alignment_locus in partial_matches and alignment.hit_def == alignment_locus+":1":

                for hsp in alignment.hsps:
                    # if alignment_locus == testlocus:
                    #     print("-1 is present",hsp.align_length)
                    # if testlocus in alignment_locus:
                    #     print("\n\n")
                    #     print(alignment.hit_def)
                    #     # print('allele size:',sizes[alignment_locus+"-1"])
                    #     print('alignment length:', hsp.align_length)
                    #     print('identities:', hsp.identities)
                    #     print('gaps:', hsp.gaps)
                    #     print('query contig', result.query)
                    #     print('query start', hsp.query_start)
                    #     print('query end: ', hsp.query_end)
                    #     print('subject start', hsp.sbjct_start)
                    #     print('subject end: ', hsp.sbjct_end)
                    #     print(hsp.query)
                    #     print(hsp.match)
                    #     print(hsp.sbjct)
                    if alignment_locus not in partials:
                        partials[alignment_locus] = [(hsp,result.query,alignment.hit_def)]
                    else:
                        partials[alignment_locus].append((hsp,result.query,alignment.hit_def))

                    # for queries with non intact hit to ref allele check other hits to see if remainder of locus is present before calling deletion - could be caused by assembly error or repetitive insertion
                    # need to deal with split matches that overlap - STMMW_45221-1 in current example 100 gene scheme
    rmlis = []

    partials2 = {}
    hitref = {}
    for locus in partials:
        hsplis = []
        c = 1
        ## if 2 or more hsps pass filter but overlap by more that 60% call as 0
        if len(partials[locus]) > 1:
            olcheck = check_for_multiple_ol_partial_hsps(partials[locus],hsp_ident)
            if olcheck == "overlap":
                rmlis.append(locus)
                # print("hits overlap by >60%",locus)
                # for hit in partials[locus]:
                #     print(hit[0],"\n\n\n")
            else:
                locuspartials = remove_hsps_entirely_within_others(partials[locus],0,hsp_ident)
        else:
            locuspartials = partials[locus]
            olcheck = "no overlap"
        if olcheck =="no overlap":
            locuspartials2 = check_ends_for_snps(locuspartials, ref_alleles[locus],"",locus,qgenome)
            partials2[locus] = locuspartials
            # if testlocus == locus:
            #     for hsp in partials2[locus]:
            #         # print(hsp)
            #         print_blasthsp("", "", hsp[0], "")
            if len(locuspartials2) == 0:
                rmlis.append(locus)
            # if len(partials[locus]) > 0 and len(locuspartials) == 0:
                # print(locus)
                # hsp = partials[locus][0][0]
                # print_blasthsp("","",hsp,"")
                # print(len(partials[locus]))
                # print(len(locuspartials))

            # for hsp in locuspartials:
            #     hs = hsp[0].query.replace("-","")
            #     hitref[locus+"_hit_no_"+str(c)] = hs
            #     hitseq = Seq(hs)
            #     # hitseq = only_nucs(hitseq)
            #     s = SeqRecord.SeqRecord(hitseq,locus+"_hit_no_"+str(c),description="")
            #     hsplis.append(s)
            #     c+=1
            # if len(locuspartials) != 0:
            #     SeqIO.write(hsplis,"tmp/"+locus+"_hits.fasta","fasta")

    npartials = {}
    for locus in partials2:
        if locus not in rmlis:
            npartials[locus] = partials2[locus]
            # if locus == "STMMW_44761":
            #     for p in npartials[locus]:
            #         print(p[0])

    # if "STMMW_40231" in npartials:
    #     print(npartials["STMMW_40231"])


    return npartials,hitref


def get_combined_hsp_coverage_of_ref_allele(reflen,hspls,hsp_thresh):
    range_lis = []
    for tup in hspls:
        hsp = tup[0]
        if hsp_filter_ok(hsp,30,hsp_thresh):
            if hsp.sbjct_start > hsp.sbjct_end:
                range_lis.append((hsp.sbjct_end, hsp.sbjct_start))
            else:
                range_lis.append((hsp.sbjct_start,hsp.sbjct_end))
    mrange_lis = utils.merge_intervals(range_lis)
    tot_covered = 0
    for hit in mrange_lis:
        tot_covered += hit[1]-hit[0]
    fraction_covered = float(tot_covered)/float(reflen)
    return fraction_covered


def check_ends_split_contigs(hsp_start_tup,hsp_end_tup,reflen,alleleseq):
    start_contig = hsp_start_tup[1]
    end_contig = hsp_end_tup[1]
    start_hsp = hsp_start_tup[0]
    end_hsp = hsp_end_tup[0]

    #start hsp

    new_qstart = 0
    new_qend = 0
    sstart = start_hsp.sbjct_start
    send = start_hsp.sbjct_end
    qstart = start_hsp.query_start
    qend = start_hsp.query_end
    added_start = ""

    if int(sstart) > int(send):
        added_start = "N" * (int(reflen) - int(sstart))
        # new_qstart = int(qstart) - (reflen - int(sstart))
        # added_st = start_contig[new_qstart:int(qstart) - 1]
        # if added_st.count("N") == len(added_st):
        #     added_start = reverse_complement(added_st)
        # else:
        #     added_start = "N"*len(added_st)
    else:
        added_start = "N"*int(sstart-1)
        # new_qstart = int(qstart) - int(sstart)
        # added_st = start_contig[new_qstart:int(qstart) - 1]
        # if added_st.count("N") == len(added_st):
        #     added_start = added_st
        # else:
        #     added_start = "N"*len(added_st)

    #end_hsp

    sstart2 = end_hsp.sbjct_start
    send2 = end_hsp.sbjct_end
    qstart2 = end_hsp.query_start
    qend2 = end_hsp.query_end

    added_end = ""

    if int(sstart2) > int(send2):
        added_end = "N" * (int(send) - 1)
        # new_qend = int(qend2) + int(send2) - 1
        # added_e = end_contig[int(qend2):new_qend - 1]
        # if added_e.count("N") == len(added_e):
        #     added_end = reverse_complement(added_e)
        # else:
        #     added_end = "N"*len(added_e)
    elif int(sstart2) < int(send2):
        added_end = "N" * (reflen - int(send))
        # new_qend = int(qend2) + (reflen - int(send2) - 1)
        # added_e = end_contig[int(qend2):new_qend - 1]
        # if added_e.count("N") == len(added_e):
        #     added_end = added_e
        # else:
        #     added_end = "N"*len(added_e)

    full_allele = added_start + alleleseq + added_end

    return full_allele, added_start, added_end


def check_split_over_contigs(hspls,query_genome,reflen):
    """
    Check if matches are at ends of contigs (or are that last non N part of the contigs)
    if so get any overlap that may have occured

    :param contigs_dict: dict of {contig_name:hsp(s) for that contig}
    :param query_genome:
    :param reflen:
    :return:
    """
    range_lis = []
    range_d = {}
    for tup in hspls:
        hsp = tup[0]
        # print(tup[1])
        # print('alignment length:', hsp.align_length)
        # print('identities:', hsp.identities)
        # print('gaps:', hsp.gaps)
        # print('ref_allele_len',str(reflen))
        # print('query start', hsp.query_start)
        # print('query end: ', hsp.query_end)
        # print('subject start', hsp.sbjct_start)
        # print('subject end: ', hsp.sbjct_end)
        # print(hsp.query)
        # print(hsp.match)
        # print(hsp.sbjct,"\n")
        if hsp.sbjct_start > hsp.sbjct_end:
            sst = hsp.sbjct_end
            send = hsp.sbjct_start
        else:
            sst = hsp.sbjct_start
            send = hsp.sbjct_end
        t = (sst,send)
        range_lis.append(t)
        range_d[t] = tup


    sorted_range_list = sorted(range_lis,key=lambda tup: tup[0])
    hsp1 = range_d[sorted_range_list[0]][0]

    newseq = hsp1.query

    if hsp1.sbjct_start > hsp1.sbjct_end:
        newseq = utils.reverse_complement(hsp1.query)

    for rang in range(len(sorted_range_list)-1):
        cur_range_end = sorted_range_list[rang][1]
        next_range_start = sorted_range_list[rang+1][0]
        cur_hsp = range_d[sorted_range_list[rang]][0]
        next_hsp = range_d[sorted_range_list[rang+1]][0]
        if cur_range_end > next_range_start-1:
            overlap = cur_range_end - next_range_start
            #last part of curhsp and 1st part of next hsp

            cur_query = cur_hsp.query
            if cur_hsp.sbjct_start > cur_hsp.sbjct_end:
                cur_query = utils.reverse_complement(cur_query)

            next_query = next_hsp.query
            if next_hsp.sbjct_start > next_hsp.sbjct_end:
                next_query = utils.reverse_complement(next_query)

            end_cur = cur_query[-1*(overlap+1):]
            start_next = next_query[:overlap+1]
            if end_cur == start_next:
                newseq += next_query[overlap+1:]
            else:
                return "inconsistent overlap",""
            #check overlap for identity
        elif cur_range_end == next_range_start-1:
            add = next_hsp.query
            if next_hsp.sbjct_start > next_hsp.sbjct_end:
                add = utils.reverse_complement(next_hsp.query)
            newseq += add
        elif cur_range_end < next_range_start-1:
            #add Ns in gap
            missing_no = next_range_start - cur_range_end-1
            nstring = "N"*missing_no
            newseq += nstring
            add_seq = next_hsp.query
            if next_hsp.sbjct_start > next_hsp.sbjct_end:
                add_seq = utils.reverse_complement(next_hsp.query)
            newseq += add_seq
        else:
            print("\n\nPROBLEM SPLIT CONTIGS\n\n")
    if len(newseq) == reflen:
        return "reconstructed split",newseq
    else:
        st_hsp = range_d[sorted_range_list[0]]
        en_hsp = range_d[sorted_range_list[-1]]
        full_allele, added_start, added_end = check_ends_split_contigs(st_hsp,en_hsp,reflen,newseq)
        # print("FULL ALLELE\n",full_allele,"\n\n")
        return "reconstructed split w ends", full_allele
        # TODONE RUN check ends - may need to modify for different hsps and contigs at ends
        # - wrote check_ends_split_contigs

    # TODONE see below
    #first sort by start number
    #if overlap check overlap region for differences
    #if diff return multiple copies
    #else
    # if start and end are complete return allele as "pseudointact"

    #if no overlap add Ns to internal


def check_all_orient(hsp_list):
    """

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


def check_matching_overlap(hsp1,hsp2,orient):
    if orient == "positive":
        ol_len = hsp1.sbjct_end - hsp2.sbjct_start
        if hsp1.query[-1*ol_len:] != hsp2.query[:ol_len]:
            return "nomatch",ol_len
        else:
            return "match",ol_len
    elif orient == "negative":
        ol_len = hsp1.sbjct_start - hsp2.sbjct_end
        if hsp1.query[:ol_len] != hsp2.query[-1*ol_len:]:
            return "nomatch",ol_len
        else:
            return "match",ol_len


def largest_nonn_strings(string):
    """

    :param string:
    :return:
    """
    matches = re.findall(r"[ATGC]*",string)
    mx = max([len(x) for x in matches])
    return mx


def check_ends(contig,qstart,qend,sstart,send,reflen,alleleseq,locus):
    """

    :param contig:
    :param qstart:
    :param qend:
    :param sstart:
    :param send:
    :param reflen:
    :param alleleseq:
    :return:
    """
    # if hit in reverse orientation have to check opposite end of hit in query genome

    new_qstart = 0
    new_qend = 0
    added_start = ''
    added_end = ''

    full_allele = ""

    if int(sstart) > int(send):
        new_qstart = int(qstart) - (reflen - int(sstart))
        new_qend = int(qend) + int(send) - 1
        added_start = "N"*(int(reflen) - int(sstart))
        #TODO Indel: following commented out elifs check for identity of missing region; if not Ns then assumed to be truncated allele - however when ignoring dels we can just add Ns to all of these regions becuase we treat dels and N regions the same
        # added_st = contig[new_qstart:int(qstart) - 1]
        # if added_st.count("N") == len(added_st):
        #     added_start = reverse_complement(added_st)
        # else:
        #     added_start = "N"*(int(reflen) - int(sstart))
        #     print("len of added start", int(reflen) - int(sstart))
        added_end = "N"*(int(send)-1)
        # added_e = contig[int(qend):new_qend - 1]
        # if added_e.count("N") == len(added_e):
        #     added_end = reverse_complement(added_e)
        # else:
        #     added_end = "N"*(int(send)-1)
        #     print("len of added end",int(send)-1)
        full_allele = added_end + alleleseq + added_start

    elif int(sstart) < int(send):
        new_qstart = int(qstart) - int(sstart)
        new_qend = int(qend) + (reflen - int(send) - 1)
        added_start = "N"*int(sstart-1)
        # added_st = contig[new_qstart:int(qstart) - 1]
        # if added_st.count("N") == len(added_st):
        #     added_start = added_st
        # else:
        #     added_start = "N"*int(sstart-1)
        added_end = "N"*(reflen - int(send))
        # added_e = contig[int(qend):new_qend - 1]
        # if added_e.count("N") == len(added_e):
        #     added_end = added_e
        # else:
        #     added_end = "N"*(reflen - int(send))
        full_allele = added_start + alleleseq + added_end
    # if locus == testlocus:
    #     print(added_start)
    #     print(added_end)


    # TODOFIXED: check that this last section actually does something
    # if missing flanking region is not all Ns then use original query hit edge - allows for rearrangements/deletions where flanking regions are really gone

    # full_allele = added_start + alleleseq + added_end

    return full_allele,added_start,added_end


def check_mid(contig,hspls,q_genome,wordsize,reflen,locus):
    """
    use word length ratio as cutoff for allowing non-N letters in restored gap sequence(i.e. if word length is 12. non-N sequence chunk can be a max of 18 before mid regions is called something else)
    :param contig:
    :param hspls:
    :param q_genome:
    :return:
    """


    orient = check_all_orient(hspls)
    # if locus == testlocus:
    #     for hsp in hspls:
    #         print(hsp)
    if orient == "mixed":
        # TODO work out what to do with mixed orientation, currently call as 0
        return "mixed_orientation"
    else:
        #TODONE work out possible issue with orientation being messed up
        ordered_by_query = []
        order = {}
        for tup in hspls:
            hsp = tup[0]
            order[hsp.query_start] = hsp
            # if locus == testlocus:
            #     print_blasthsp("","",hsp,"")

        sorted_hsps = sorted(map(int, order.keys()))
        # hsp1 = order[sorted_hsps[0]]


        full_allele = remove_indels_from_hsp(order[sorted_hsps[0]])#order[sorted_hsps[0]].query
        sstart = order[sorted_hsps[0]].sbjct_start
        send = order[sorted_hsps[-1]].sbjct_end
        qstart = order[sorted_hsps[0]].query_start
        qend = order[sorted_hsps[-1]].query_start

        if orient == "positive":
            for start in range(len(sorted_hsps)-1):
                hsp = order[sorted_hsps[start]]
                hspnext = order[sorted_hsps[start+1]]
                mid_section_st = int(hsp.query_end)
                mid_section_en = int(hspnext.query_start)
                if mid_section_en < mid_section_st:
                    olcheck, ollen = check_matching_overlap(hsp, hspnext, "negative")
                    add = q_genome[contig][hspnext.query_start+ollen:hspnext.query_end-1]
                    full_allele = full_allele + add
                    # mid_section_seq = "N"*(hspnext.sbjct_start-hsp.sbjct_end-1)
                else:
                    mid_section_seq = "N"*(hspnext.sbjct_start-hsp.sbjct_end-1)
                    # mid_section_seq = q_genome[contig][mid_section_st:mid_section_en-1]
                    # if locus == testlocus:
                    #     print(order[sorted_hsps[start+1]])
                    full_allele += mid_section_seq
                    full_allele += remove_indels_from_hsp(hspnext)#hspnext.query
                    nonn_size = largest_nonn_strings(mid_section_seq)
                    if nonn_size > (wordsize*2):
                        return "possible_insertion"
        elif orient == "negative":
            for start in range(len(sorted_hsps)-1):
                hsp = order[sorted_hsps[start]]
                hspnext = order[sorted_hsps[start+1]]
                mid_section_st = int(hsp.query_end)
                mid_section_en = int(hspnext.query_start)
                if mid_section_en < mid_section_st:
                    olcheck,ollen =check_matching_overlap(hsp,hspnext,"negative")
                    add = q_genome[contig][hspnext.query_start+ollen:hspnext.query_end-1]
                    full_allele = add + full_allele
                    # mid_section_seq = "N" * (hsp.sbjct_end - hspnext.sbjct_start-1)
                else:
                    mid_section_seq = "N"*(hsp.sbjct_end - hspnext.sbjct_start-1)
                    # mid_section_seq = q_genome[contig][mid_section_st:mid_section_en - 1]
                    full_allele = mid_section_seq + full_allele
                    full_allele = remove_indels_from_hsp(hspnext) + full_allele#hspnext.query
                    nonn_size = largest_nonn_strings(mid_section_seq)
                    if nonn_size > (wordsize*2):
                        return "possible_insertion"

        full_allele, added_start, added_end = check_ends(contig, qstart, qend, sstart, send, reflen,full_allele,locus)
        # if locus == testlocus:
        #     print("fixmid")
        #     print("start " + added_start)
        #     print("end " + added_end)
        #     print("mid " + mid_section_seq)
        #     print("full " + full_allele)
        return full_allele


def generate_query_allele_seqs(partial_hsps,query_genome,alleles_sizes,missing_perc_cutoff,wordsize,ref_alleles,qgenome,hsp_thresh):
    query_genome = SeqIO.parse(query_genome,"fasta")
    q_genome = {}
    uncallable = {}
    calls = {}

    full_allele = ""

    for s in query_genome:
        q_genome[s.id] = str(s.seq)

    #check that number of identities in blast hits is at least X fraction of normal reference allele length
    # need to check that at least x of allele is covered
    #TODONE currently only counts ident numbers, could count idents to same region 50 times and come up with "intact" allele - DONE in get_combined_hsp_coverage_of_ref_allele function
    for locus in partial_hsps:
        hspls = partial_hsps[locus]
        reflen = alleles_sizes[locus+":1"]
        frac_covered = get_combined_hsp_coverage_of_ref_allele(reflen,hspls,hsp_thresh)
        if frac_covered < float(missing_perc_cutoff):
            # print(locus," unscoreable: ",frac_covered)
            uncallable[locus] = "unscorable_too_much_missing"
        # if locus == testlocus:
        #     print(frac_covered,reflen)
        # TODONE add function to remove hits that are within others i.e. best hit  = 1-1456, next hit = 120-200 -> remove 120-200 hit
        # Also removes hsps <95% identity to ref and <30bp
        hspls = list(remove_hsps_entirely_within_others(hspls,locus,hsp_thresh))
        hspls = check_ends_for_snps(hspls,ref_alleles[locus],"","",qgenome)


        hspls2 = []
        if len(hspls) == 0:
            uncallable[locus] = "failed_filter"
        elif len(hspls) > 1:
            # if locus == testlocus:
            #     print("more than one")
            contigs = {}
            hsp_to_investigate = []
            for tup in hspls:
                hsp = tup[0]
                if hsp_filter_ok(hsp,30,hsp_thresh):
                    hspls2.append(tup)
                    contig = tup[1].split(" ")[0]
                    if contig not in contigs:
                        contigs[contig] = [tup]
                    else:
                        contigs[contig].append(tup)

            # if hits come from > 1 contig
            if len(contigs.keys()) > 1:

                message,full_allele = check_split_over_contigs(hspls2,query_genome,reflen)
                # if locus == testlocus:
                #     print(message,full_allele)
                if message == "inconsistent overlap":
                    uncallable[locus] = message
                else:
                    calls[locus] = full_allele
                # TODONE Write split over contigs function
                # return message,full_allele
            else:
                # if locus == testlocus:
                #     print("one contig")
                for c in contigs:
                    fixed_mid = check_mid(c, contigs[c],q_genome,wordsize,reflen,locus)
                    # if locus == testlocus:
                    #     print(fixed_mid)
                    #     sl(1)
                    if fixed_mid == "possible_insertion":
                        uncallable[locus] = "possible_insertion"
                    elif fixed_mid == "mixed_orientation":
                        uncallable[locus] = "mixed_orientation"
                    else:
                        calls[locus] = fixed_mid
                        # if locus == "STMMW_18251" or locus == "STM3980":
                        #     print("mid",full_allele)
                        # print("mid")

        else:
            # if locus == testlocus:
            #     print("one hit only")
            #     sl(1)
            # print(locus+": "+hspls[0][1])
            hsp = hspls[0][0]
            queryid = hspls[0][1].split(" ")[0]
            contig = q_genome[queryid]
            # print("\n")
            # print('allele length:', reflen)
            # print('alignment length:', hsp.align_length)
            # print('identities:', hsp.identities)
            # print('gaps:', hsp.gaps)
            # print('query contig', queryid)
            # print('query start', hsp.query_start)
            # print('query end: ', hsp.query_end)
            # print('subject start', hsp.sbjct_start)
            # print('subject end: ', hsp.sbjct_end)
            # print(hsp.query)
            # print(hsp.match)
            # print(hsp.sbjct)
            seq = remove_indels_from_hsp(hsp)
            # print(seq)
            full_allele,added_start,added_end = check_ends(contig,hsp.query_start,hsp.query_end,hsp.sbjct_start,hsp.sbjct_end,reflen,seq,locus)
            # if locus == testlocus:
            #     print("one hit only")
            #     print(queryid)
            #     print("start "+added_start)
            #     print("end "+added_end)
            #     print("full "+full_allele)
            #     sl(5)
            calls[locus]=full_allele
    calls2 = dict(calls)
    for locus in calls:
        if float(len(calls[locus])) > 1.5*float(alleles_sizes[locus+":1"]):
        # print(locus,"unscorable: N-locus too long!!")
            uncallable[locus] = "unscorable_too_long"
        else:
            calls2[locus] = calls[locus]
            # if locus == testlocus:
            #     print(calls2[locus])

    # if "STMMW_40231" in calls2:
    #     print(calls2["STMMW_40231"])

    return calls2,uncallable

def invert_dict(d):
    inverse = dict()
    for key in d:
        # Go through the list that is saved in the dict:
        for item in d[key]:
            # Check if in the inverted dict the key exists
            if item not in inverse:
                # If not create a new list
                inverse[item] = [key]
            else:
                inverse[item].append(key)
    return inverse


def mask_high_snp_regions(locus,recon_locus,ref_locus,window_size,snp_limit):

    if len(recon_locus) != len(ref_locus):
        return recon_locus

    halfwindow = int(window_size/2)
    mutpos = []
    outlocus = list(str(recon_locus))
    for pos in range(len(recon_locus)):
        if recon_locus[pos] != ref_locus[pos] and recon_locus[pos] not in ["N", "-"] and ref_locus[pos] not in ["N", "-"]:
            mutpos.append("X")
        else:
            mutpos.append("M")

    for x in range(halfwindow,len(ref_locus)-halfwindow):
        window = mutpos[x-halfwindow:x+halfwindow]
        if window.count("X") > snp_limit:
            for pos in range(x-halfwindow,x+halfwindow):
                outlocus[pos] = "N"
    outlocus = "".join(outlocus)
    # if outlocus != recon_locus:
    #     print("\n\n\n")
    #     print(locus)
    #     print("\n")
    #     print(len(recon_locus),recon_locus)
    #     print("\n")
    #     print(len(mutpos),"".join(mutpos))
    #     print("\n")
    #     print(len(outlocus),outlocus)
    #     print("\n\n\n")
    return outlocus


def check_locus_allele_freqs_in_higher_cc(st_to_cc_unmerge,profdict,prev_cc_unmerge,locus_profile_position,locus,hgt_dict,scheme_pos,scheme,prevscheme,newnegs,orig_ccmerges):
    cc_merges = orig_ccmerges[scheme]
    st_to_cc = st_to_cc_unmerge[scheme]

    nst_to_cc = {}

    for ost in st_to_cc:
        occ = st_to_cc[ost]
        if occ in cc_merges:
            nst_to_cc[ost] = cc_merges[occ]
        else:
            nst_to_cc[ost] = occ
    # print(scheme_pos)
    # print(scheme)
    if scheme_pos == 0:
        return "0"
    else:
        prevmerges = orig_ccmerges["MGT"+str(scheme_pos+1)]
        # print("MGT"+str(scheme_pos+1))
        if prev_cc_unmerge in prevmerges:
            prev_cc = prevmerges[prev_cc_unmerge]
        else:
            prev_cc = prev_cc_unmerge


    locus_index = locus_profile_position[scheme][locus]
    cc_to_st = dict(invert_dict(nst_to_cc))

    #get clonal complex of higher scheme (n-1) in current strain then get all CC in scheme n that have that n-1 scheme
    other_cc_list = []

    for i in hgt_dict:
        hgt = i.split("-")
        # print(hgt)
        if hgt[prevscheme] == prev_cc:
            other_cc_list.append(hgt[scheme_pos])

    #get all sts that match those scheme n CCs
    stlis = []
    # print(other_cc_list)
    for cc in other_cc_list:
        if cc in cc_to_st:
            stlis += cc_to_st[cc]
    # stlis = [tuple([x[0],x[1]x.split(".")) for x in stlis]
    # print(stlis)
    #using locus of interest index get allele ids for above sts
    allele_lis = []
    for i in profdict[scheme]:
        if i in stlis:
            # print(profdict[scheme][i])
            # print(len(profdict[scheme][i]),locus_index)
            allele_lis.append(profdict[scheme][i][locus_index])

    #get allele ids that match current new allele
    called_alleles = [x[1] for x in newnegs]

    #iterate over new allele matches and related strain allele matches, look for any match
    match = []
    for i in allele_lis:
        if i != "0":
            stold = i.split("_")[0].strip("-")
            for j in called_alleles:
                stnew = j.split("_")[0].strip("-")
                if stold == stnew:
                    match.append(stold)

    match = [x for x in match if x != "0"]
    # print(allele_lis)
    # print(called_alleles)
    # print(match)
    #if no matches return 0
    if match == []:
        return "0"
    #if one st match (any count) return that
    elif len(list(set(match))) == 1:
        # print("negs",match[0])
        return match[0]
    #if more than one match return most common
    else:
        outallele = utils.most_common(match)
        # print("negs", outallele)
        return outallele


def get_neg_matches(existing,reconstructed,sizes,Nfraction,st_to_cc,profdict,prev_cc,locus_profile_position,hgt_dict,scheme_pos,scheme,prevscheme,orig_ccmerges,conf):
    #loci with perfect matches to intact alleles and to reference are already removed
    #reconstructed alleles have Ns where indels are (TODOindel reconstructed will need to have "-" and "N" to distinguish masked and indel)
    # print("pre", (" --- %s seconds ---" % (time.time() - start_time)))
    outcomes = {}
    # print(len(reconstructed))

    snp_density_win = int(conf["MGT_conf"]["snp_density_filter_window"])
    snp_density_max_num = int(conf["MGT_conf"]["snp_density_filter_limit"])

    for newlocus in reconstructed:
        newlocus = str(newlocus)
        outcomes[newlocus] = ()
        for existlocus in existing[newlocus]:
            if reconstructed[newlocus] == existing[newlocus][existlocus]:
                ## at this point can get reconstructed alleles that match intact pre-existing
                if "N" in reconstructed[newlocus]:
                    outcomes[newlocus] = ("negmatch", existlocus, "",[])
                else:
                    outcomes[newlocus] = ("posmatch",existlocus,"",[])

                # print("CHECK1")
                # print(">1\n"+reconstructed[newlocus])
                # print(">2\n"+existing[newlocus][existlocus])
                # print("CHECK2\n")
        #locis with exact match to existing neg or pos are removed above
        #only loci with no exact matches remain - can be new pos, novel neg, or new existing neg
        muts = []
        newnegs = []
        if outcomes[newlocus] == ():
            newseq = str(reconstructed[newlocus])
            # if newlocus == testlocus:
            #     print(newlocus)
            #     print(newseq)
            #     print(existing[newlocus]["1"])
            #     print("\n")
            newseq = mask_high_snp_regions(newlocus,newseq,existing[newlocus]["1"],snp_density_win,snp_density_max_num)

            if len(newseq) != sizes[newlocus+":1"]:
                outcomes[str(newlocus)] = ("new pos allele", "0", newseq, muts)
                # print(newlocus, len(newseq), sizes[newlocus+":1"])
                # print(newseq)
                # print(existing[newlocus]["1"])
            elif newseq.count("N") > Nfraction*sizes[newlocus + ":1"]:
                outcomes[newlocus] = ("new pos allele", "0", newseq, muts)
            else:
                # print(newlocus)
                oldseq = existing[newlocus]["1"]
                for pos in range(sizes[newlocus + ":1"]):
                    if newseq[pos] != oldseq[pos] and newseq[pos] not in ["N", "-"] and oldseq[pos] not in ["N", "-"]:
                        muts.append((oldseq[pos], str(pos), newseq[pos]))

                for existlocus in existing[newlocus]:
                    oldseq = str(existing[newlocus][existlocus])
                    anymut = 'no'
                    for pos in range(sizes[newlocus+":1"]):
                        if newseq[pos] != oldseq[pos] and newseq[pos] not in ["N","-"] and oldseq[pos] not in ["N","-"]:
                            anymut = 'yes'
                            # if newlocus == testlocus:
                            #     print("!!!!!!!!new",newseq[pos],"old",oldseq[pos], "exist_allele: ",existlocus)
                    if anymut == 'no':
                        # if newlocus == testlocus:
                        #     print("no snps: ",existlocus)
                        #any loci at this point must be neg matches if no muts
                        allele = existlocus.split("_")[0].strip("-")
                        newnegs.append(("new neg allele", allele, newseq,muts))

                ##need to check if neg match is contradicted by pos match and remove neg match if so.

                # print(newnegs)
                posmatch = []
                matchallelelist = [x[1] for x in posmatch]
                nnewnegs = []
                for match in newnegs:
                    if "-" in match[1]:
                        posversion = match[1].split("_").strip("-")
                        if posversion in existing[newlocus]:
                            if posversion not in matchallelelist:
                                continue
                            else:
                                nnewnegs.append(match)
                        else:
                            nnewnegs.append(match)
                    else:
                        nnewnegs.append(match)

                newnegs = list(nnewnegs)

                # print(newnegs)
                # print("\n")

                if len(newnegs) == 0:
                    if "N" in newseq:
                        # print("0hit", newlocus, newnegs)
                        ##call as 0
                        outcomes[newlocus] = ("novel neg allele", "0", newseq, muts)
                    else:
                        outcomes[newlocus] = ("novel pos allele", "", newseq, muts)
                elif len(newnegs) == 1:
                    # print("1hit",newlocus,newnegs)
                    if "N" in newseq:
                        allele = check_locus_allele_freqs_in_higher_cc(st_to_cc,profdict,prev_cc,locus_profile_position,newlocus,hgt_dict,scheme_pos,scheme,prevscheme,newnegs,orig_ccmerges)
                        ##check what other related sts have at this position (by higher CC) if any called with one that can match then assign otherwise 0
                        outcomes[newlocus] = ("new neg allele", allele, newseq, muts)
                    else:
                        ## check if negative matches have any pos allele associated if no then assign to neg number - done later
                        allele = newnegs[0][1]
                        outcomes[newlocus] = ("new pos allele", allele, newseq, muts)
                else:
                    allele_hits = list(set([x[1] for x in newnegs]))
                    # print(allele_hits)
                    # print(">2hit",newlocus,newnegs)
                    if len(allele_hits) > 1:
                        if "N" in newseq:
                            allele = check_locus_allele_freqs_in_higher_cc(st_to_cc, profdict, prev_cc, locus_profile_position,
                                                                  newlocus, hgt_dict, scheme_pos, scheme,prevscheme,newnegs,orig_ccmerges)
                            ##check what other related sts have at this position (by higher CC) if any called with one that can match then assign otherwise 0
                            outcomes[newlocus] = ("new neg allele", allele, newseq, muts)
                        else:
                            ## check if negative matches have any pos allele associated if no then assign to neg number otherwise next pos- done later
                            outcomes[newlocus] = ("new pos allele", "", newseq, muts)
                    else:
                        if "N" in newseq:
                            allele = check_locus_allele_freqs_in_higher_cc(st_to_cc, profdict, prev_cc, locus_profile_position,
                                                                  newlocus, hgt_dict, scheme_pos, scheme,prevscheme,newnegs,orig_ccmerges)
                            ##check what other related sts have at this position (by higher CC) if any called with one that can match then assign otherwise 0
                            outcomes[newlocus] = ("new neg allele", allele, newseq, muts)
                        else:
                            ## check if negative matches have any pos allele associated if no then assign to neg number that is most frequent in clade - done later
                            allele = newnegs[0][1]
                            outcomes[newlocus] = ("new pos allele", allele, newseq, muts)
        # print(newlocus,outcomes[newlocus][0],outcomes[newlocus][1],outcomes[newlocus][3])
        # print(reconstructed[newlocus])
        # print("\n")
        # print(existing[newlocus]["1"],"\n\n")
        # sl(1)
    # print("post", (" --- %s seconds ---" % (time.time() - start_time)))
    # for i in outcomes:
    #     if outcomes[i][0] == "posmatch" and "N" in str(reconstructed[i]):
    #         print("ERROR for "+ i)
    # if "STMMW_40231" in outcomes:
    #     print(outcomes["STMMW_40231"])

    return outcomes


def assign_new_allele_names(locus,type,negallele,existing_alleles):
    """

    :param locus: locus being named
    :param type: positive, negative or novel negative - what type of assignment is happening
    :param negallele: if it is a neg allele which pos allele is it based on
    :param existing_alleles: list of existing alleles for that locus
    :return: next allele number in the sequence
    """
    posalleles = [x.split("_")[0].strip("-") for x in existing_alleles ]
    posalleles = [x for x in posalleles if x.isdigit()]



    if type == "novel pos allele":

        newallele = max(map(int,posalleles))+1

        return newallele
        # assign new pos allele i.e. 8 where 7 is current highest

    elif type ==  "new pos allele":
        newposallele = negallele
        if newposallele not in posalleles:
            newallele = newposallele
        else:
            newallele = max(map(int, posalleles)) + 1
        return newallele

    elif type == "new neg allele":
        neglis = [x for x in existing_alleles if "-" in x]
        neglis = [x for x in neglis if x.split("_")[0][1:] == negallele]
        neglis = [x.split("_")[-1] for x in neglis]
        # neglis = [x.split("_")[-1] for x in existing_alleles if x.split(":")[0].strip("-") == negallele]
        # if locus == testlocus:
        #     print(existing_alleles)
        #     print(neglis)

        if len(neglis) > 0:

            newnegnumber = max(map(int, neglis)) + 1
            newallele = "-"+negallele+"_"+str(newnegnumber)

            return newallele

        else:
            newallele = "-" + negallele + "_1"
            return newallele
        # assign new neg i.e. -2_10 where -2_9 was highest -2 allele

    elif type == "novel neg allele":

        newallele = max(map(int, posalleles)) + 1
        newallele = "-"+str(newallele)+"_1"

        return newallele
        # assign neg allele to new overall allele i.e. -8.1 where 8 does not exist


def call_alleles(schemelist, alleles, query_genome, locus_refrence_locs, scheme, wordsize, locuscalls, all_alleles,
                 subset, indels, ref_called_alleles, ref_alleles_fastas, ref_blast_hits,conf,qgenome,st_to_cc,profdict,
                 prev_cc,locus_profile_position,hgt_dict,global_schemels,prevscheme,orig_ccmerges):

    missing_limit = float(conf["MGT_conf"]["min_allele_not_missing"])
    print(scheme)
    start_time1 = time.time()

    """

    :param schemelist:
    :param alleles:
    :param query_genome:
    :param locus_refrence_locs:
    :param scheme:
    :param wordsize:
    :param locuscalls:
    :param all_alleles:
    :param subset:
    :return:
    """

    if os.path.exists("tmp"):
        shutil.rmtree('tmp')
        os.mkdir("tmp")
    else:
        os.mkdir("tmp")

    hsp_ident_thresh = float(conf["MGT_conf"]["hsp_identity_threshold"])

    alleles_folder = alleles + "/*.fasta"

    # gets list of allele fasta file paths
    allele_fastas = IO.get_allele_fastas(alleles_folder)

    new_allele_outdict = {}

    # gets list of loci in scheme of interest

    schemelist1 = list(schemelist)
    schemelist2 = list(schemelist)
    schemelist3 = list(schemelist)
    schemelist4 = list(schemelist)

    scheme_fastas = []

    # makes a list of all of the fasta files in the scheme of interest
    # also makes list of word size to use depending on ref allele size
    c = 0
    c1 = 0

    for i in schemelist2:
        c += 1
        if i not in ref_called_alleles:
            allele_file = allele_fastas[i]
            scheme_fastas.append(allele_file)
            # ws = get_word_size(allele_file)
            # wordsizes.append(ws)
            c1 += 1
    print("p1 %s " % (time.time() - start_time1))
    start_time1 = time.time()
    if len(scheme_fastas) == 0:
        for i in schemelist2:
            locuscalls[i] = ref_called_alleles[i]
        return locuscalls, new_allele_outdict, indels

    # writes all fasta files in scheme to tmpdb location -- effectively a python "cat"
    tmpdb = "tmp/" + scheme + "_alleles.fasta"


    existing_alleles, ref_alleles, allele_seqs = make_intact_allele_fasta(scheme_fastas, tmpdb, all_alleles, subset,
                                                                          scheme)
    print("p2 - intact_allele_fasta %s" % (time.time() - start_time1))
    start_time1 = time.time()

    # load scheme allele fastas as Seq object
    allele_sizes = {}
    allele_seqs_byloc = {}


    # make dict of allele sizes
    for allele in allele_seqs:
        allele_sizes[allele] = len(allele_seqs[allele])
        loc = str(allele.split(":")[0])
        all = str(allele.split(":")[1])

        if loc not in allele_seqs_byloc:
            allele_seqs_byloc[loc] = {all: allele_seqs[allele]}
        else:
            allele_seqs_byloc[loc][all] = allele_seqs[allele]

    hit_locs = IO.make_contig_dict(query_genome)

    # print("preblast", (" --- %s seconds ---" % (time.time() - start_time)))
    blast_hits = wrappers.run_blast(query_genome, tmpdb, wordsize, 10, 90)
    # print("postblast", (" --- %s seconds ---" % (time.time() - start_time)))
    print("p3 - runblast %s" % (time.time() - start_time1))
    start_time1 = time.time()
    # complete loci are exact matches to an allele already in the db
    # partial loci are the rest
    partial_loci, complete_loci, hit_locs, indels = get_exactmatches(blast_hits, schemelist1, allele_sizes, hit_locs,
                                                                     indels, ref_called_alleles,hsp_ident_thresh)


    print("p3 - get_exactmatches %s" % (time.time() - start_time1))
    start_time1 = time.time()
    # missing_limit = 0.5
    # returns hsps that give the positions of matches in the query genomes that match the "1"/reference allele
    # also writes these regions to one fasta file for each locus different hsps are called "_hit_no_1", "hit_no_2" etc
    partial_hsps, hitref = get_partial_match_query_region(ref_blast_hits, partial_loci, allele_sizes, ref_alleles,qgenome,hsp_ident_thresh)

    print("p4 - get_partial_match_query_region %s" % (time.time() - start_time1))
    start_time1 = time.time()

    reconstructed, uncallable = generate_query_allele_seqs(partial_hsps, query_genome, allele_sizes, missing_limit,
                                                           wordsize, ref_alleles,qgenome,hsp_ident_thresh)

    print("p5 - generate_query_allele_seqs %s" % (time.time() - start_time1))
    start_time1 = time.time()

    scheme_pos = global_schemels.index(scheme)

    outcome = get_neg_matches(allele_seqs_byloc, reconstructed, allele_sizes, missing_limit,st_to_cc,profdict,prev_cc,locus_profile_position,hgt_dict,scheme_pos,scheme,prevscheme,orig_ccmerges,conf)

    print("p6 - get_neg_matches %s" % (time.time() - start_time1))
    start_time1 = time.time()
    # run blast for each locus with 'query genome region that corresponds to "reference" allele hit regions' against all alleles for that locus

    # todone insert function here to make intact alleles and use them to blast allele db
    # print("pre2blast", (" --- %s seconds ---" % (time.time() - start_time)))
    # neg_matches,sec_results,refhits = exactmatch_or_sec_blast(reconstructed,allele_fastas,ref_alleles_fastas)
    # print("post2blast", (" --- %s seconds ---" % (time.time() - start_time)))

    '''sec_results = run_secondary_n_blast(partial_hsps,allele_fastas)
'''
    # -----------------
    # different process for split hits
    # check for Ns
    # check for missing ends
    # check for indels - always call relative to existing allele
    #       if gap in subject = insertion in query
    #       if gap in query = deletion in query
    # check for SNPs
    # ------------------

    # if match against reference runs from reference start to end then snps/indesl can be called without checking for missing ends
    # therefore can extract snps from alignment against "-1" reference allele using locus blast results - find_snps function

    # print(scheme,schemelist)
    unacounted = list(schemelist4)
    # print(len(unacounted))
    allele_assignments = {}

    # print("exact match")

    for locus in complete_loci:
        # if locus == "STMMW_29661":
        #     print("complete",complete_loci[locus])
        matchname = complete_loci[locus]
        matchallele = matchname.split(":")[1]
        # print(matchname,matchallele)
        allele_assignments[locus] = str(matchallele)  # (matchallele,"existing")
        unacounted.remove(locus)
        # print(locus,matchallele)

    # ("new_exist_neg", allele, newseq,muts)

    new_alleles = {}
    new_allele_outdict = {}
    all_muts = {}
    for locus in outcome:
        # if locus == "STMMW_40231":
        #     print(outcome[locus])
        if outcome[locus][0] == "posmatch" or outcome[locus][0] == "negmatch":
            allele_assignments[locus] = outcome[locus][1]
            unacounted.remove(locus)
        elif outcome[locus][1] == "0":
            allele_assignments[locus] = "0"
        else:
            newno = assign_new_allele_names(locus, outcome[locus][0], outcome[locus][1], existing_alleles[locus][1])
            unacounted.remove(locus)
            allele_assignments[locus] = str(newno)
            new_alleles[locus] = outcome[locus][2]
            all_muts[locus] = outcome[locus][3]

    for locus in new_alleles:

        if locus not in uncallable:
            new_allele_outdict[locus] = (locus, allele_assignments[locus], new_alleles[locus], all_muts[locus], alleles)



    for locus in uncallable:
        allele_assignments[locus] = "0"
        if locus in unacounted:
            unacounted.remove(locus)

    # dict of allele calls {locus:allele no}
    for locus in unacounted:
        allele_assignments[locus] = "0"

    c = 0

    for locus in allele_assignments:
        if allele_assignments[locus] == "0":
            # print(locus, "ZERO")
            c += 1


    for locus in allele_assignments:
        # print(locus,allele_assignments[locus])#[0],allele_assignments[locus][1])
        # if "N" in allele_assignments[locus][2]:# and allele_assignments[locus][0] == 'new pos allele':
        #     print(locus,allele_assignments[locus])
        locuscalls[locus] = allele_assignments[locus]

    print("p7 - final_processing and output %s" % (time.time() - start_time1))

    '''
    blast hits structure: 
    list of results

    result attributes: >>>'alignments'<<<, 'application', 'blast_cutoff', 'database', 'database_length', 'database_letters', 'database_name', 'database_sequences', 'date', 'descriptions', 'dropoff_1st_pass', 'effective_database_length', 'effective_hsp_length', 'effective_query_length', 'effective_search_space', 'effective_search_space_used', 'expect', 'filter', 'frameshift', 'gap_penalties', 'gap_trigger', 'gap_x_dropoff', 'gap_x_dropoff_final', 'gapped', 'hsps_gapped', 'hsps_no_gap', 'hsps_prelim_gapped', 'hsps_prelim_gapped_attemped', 'ka_params', 'ka_params_gap', 'matrix', 'multiple_alignment', 'num_good_extends', 'num_hits', 'num_letters_in_database', 'num_seqs_better_e', 'num_sequences', 'num_sequences_in_database', 'posted_date', 'query', 'query_id', 'query_length', 'query_letters', 'reference', 'sc_match', 'sc_mismatch', 'threshold', 'version', 'window_size']
        alignment attributes: 'accession', 'hit_def', 'hit_id', >>>'hsps'<<<, 'length', 'title']
            hsp attributes: 'align_length', 'bits', 'expect', 'frame', 'gaps', 'identities', 'match', 'num_alignments', 'positives', 'query', 'query_end', 'query_start', 'sbjct', 'sbjct_end', 'sbjct_start', 'score', 'strand'

    '''

    return locuscalls, new_allele_outdict, indels