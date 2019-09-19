from time import sleep as sl
import argparse
import psycopg2
from Bio import SeqIO
import sys
from os import path

sys.path.append(path.dirname(path.dirname(path.dirname(path.abspath(__file__)))))

# from Mgt import settings

"""
make dicts with following structure:
A - {locus:{position:{mutinfo:[APidlist]}}}
B - {APid:{locus:allele}}
C - {APid:[listof strains]}
D - {locus:[list of alleles]}
to do this:
 for each isolate get MGT9 APid - make C
 for each AP for each locus get lists of alleles present - Make B
 use B to get all mut positions in a locus (from allele nos)
 for each position then iterate over each allele and assign nucleotide
 add nucleotide assignments together to make alignment


"""


def make_dash_nodash(conn, args):
    hgtquery = """ SELECT DISTINCT "locus_id" FROM "{}_allele"; """.format(args.appname)

    res = sqlquery_to_outls(conn, hgtquery)

    locuslist = [x[0] for x in res]

    dash_nodash = {x: x.replace("_", "") for x in locuslist}
    nodash_dash = {x.replace("_", ""): x for x in locuslist}

    return dash_nodash, nodash_dash


def sqlquery_to_outls(con, query):
    """
    Run sql query string with supplied connection and return results as a list of tuples
    :param con: psycopg2 sql connection object
    :param query: sql query string
    :return: sql query results (which are in a list of tuples)
    one tuple for each row returned and each tuple with the selected columns
    """
    cur = con.cursor()  # generate cursor with connection
    cur.execute(query)  # execute the sql query
    res = cur.fetchall()  # get results and save to res as list
    cur.close()  # close cursor
    return res


def get_table_nos(conn, args):
    """
    retrieve information about allele profile tables for each scheme
    :param conn: sql connection
    :param args: inputs
    :return: allele profile tables in format: {level:{table number: [list of locus names in table]}}
    """
    sqlquery = """ SELECT "scheme_id","table_num" FROM "{}_tables_ap";  """.format(args.appname)

    res = sqlquery_to_outls(conn, sqlquery)
    tables = {}
    for i in res:
        lev = int(i[0].replace("MGT", ""))

        if lev not in tables:
            tables[lev] = {i[1]: []}
        else:
            tables[lev][i[1]] = []

    for lev in tables:
        for no in tables[lev]:

            sqlcommand = "select column_name from INFORMATION_SCHEMA.COLUMNS where table_name = '{}_ap{}_{}';".format(args.appname, lev, no)

            cur = conn.cursor()

            cur.execute(sqlcommand)

            res = cur.fetchall()

            cur.close()
            # print(args.appname, lev, no)


            if lev == 9:
                if no == 0:
                    locus_columns = [x[0] for x in res[3:-6]]
                else:
                    locus_columns = [x[0] for x in res[1:-2]]
            else:
                if no == 0:
                    locus_columns = [x[0] for x in res[3:-3]]
                else:
                    locus_columns = [x[0] for x in res[1:-2]]

            tables[lev][no] = locus_columns

    return tables


def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--isolatelist", help="one per line list of isolates")
    # /Users/michaelpayne/Desktop/aus_snp_align_test_isolates.txt
    parser.add_argument("-d", "--appname", help="App name", default="Salmonella")
    parser.add_argument("-o", "--outfile", help="output fasta file")
    parser.add_argument("-f", "--fasta_output",
                        help="If used output will be a SNP alignment fasta file, if not output will be a position by position matrix",
                        action='store_true')

    args = parser.parse_args()
    return args


def get_conn():
    # database = settings.APPS_DATABASE_MAPPING[args.appname]
    #
    # args.database = database
    #
    # psql_details = settings.DATABASES[database]
    #
    # args.psqldb = psql_details['NAME']
    #
    # DbConString = "dbname='{0}' host='{1}' port='{2}' user='{3}' password='{4}'".format(psql_details['NAME'],psql_details['HOST'],psql_details['PORT'],psql_details['USER'],psql_details['PASSWORD'])  ## connection info for Db - assign new user that can only do what is needed for script
    #
    #
    # conn = psycopg2.connect(DbConString)
    # conn.autocommit = True

    DbConString = "dbname='vcseventh_15' host='0.0.0.0' port='5432' user='postgres' password='5678'"
    conn = psycopg2.connect(DbConString)
    conn.autocommit = True
    return conn


def get_apids(args, conn, maxscheme):
    missing_strains = []

    inplist = open(args.isolatelist, "r").read().splitlines()
    print(inplist)

    sqlq = """
SELECT "{appname}_isolate".identifier,"{appname}_isolate".id,"{appname}_view_apcc".ap{maxscheme}_0 
FROM "{appname}_isolate" JOIN "{appname}_view_apcc"
ON "{appname}_isolate".mgt_id = "{appname}_view_apcc".mgt_id
WHERE "{appname}_isolate".identifier in ('{idlist}')
""".format(maxscheme=maxscheme, appname=args.appname, idlist="','".join(inplist))

    res = sqlquery_to_outls(conn, sqlq)

    apidToisolate = {}
    for r in res:
        apid = r[2]
        strainid = r[1]

        print(apid, strainid)

        if apid == None:
            missing_strains.append(strainid)
        else:
            if apid not in apidToisolate:
                apidToisolate[apid] = [r[0]]
            else:
                apidToisolate[apid].append(r[0])

    return apidToisolate, missing_strains


def get_max_scheme(connection, args):
    """
    :param connection: sql connection
    :param args: inputs
    :return: highest MGT level in current database
    """
    sqlquery = """ SELECT MAX("display_order") FROM "{}_tables_ap";  """.format(args.appname)

    res = sqlquery_to_outls(connection, sqlquery)

    maxno = int(res[0][0])

    return maxno


def get_ref_pos(loc, pos, posdict):
    info = posdict[loc]
    orient = info[2]
    start = int(info[0])
    pos = int(pos)
    end = info[1]
    if orient == "+":
        refpos = start + pos
    elif orient == "-":
        refpos = end - pos
    else:
        refpos = start + pos
    return refpos


def get_apid_to_alleles(args, conn, apidToisolate, lev, table_nos, nodash_dash):
    apid_lis = list(map(str, apidToisolate.keys()))
    print(apid_lis)
    apLocAll = {}
    locus2allelels = {}
    for num in table_nos[lev]:
        tablesloci = table_nos[lev][num]

        for loc in tablesloci:
            loc = str(nodash_dash[loc])
            locus2allelels[loc] = []

        if num == 0:
            ident = "id"
        else:
            ident = "main_id"

        locuslist_string = '"' + ident + '","' + '","'.join(tablesloci) + '"'

        sqlcomm = """SELECT {0} FROM "{1}_ap{2}_{3}" WHERE {4} in ('{5}')""".format(locuslist_string, args.appname, lev,
                                                                                    num, ident, "','".join(apid_lis))

        res = sqlquery_to_outls(conn, sqlcomm)

        for i in res:
            apid = i[0]
            if apid not in apLocAll:
                apLocAll[apid] = {}
            # print(apid)
            for p in range(1, len(i)):
                # print("locus",tablesloci[p-1],"allele",i[p])
                locus = tablesloci[p - 1]
                locus = nodash_dash[locus]
                allele = i[p]
                if allele not in locus2allelels[locus]:
                    locus2allelels[locus].append(allele)
                apLocAll[apid][locus] = allele
    return apLocAll, locus2allelels


def get_snpinfo(args, conn, loc2allls):
    """
    for each allele get snp info
    for each neg allele get allele seq
    """
    # locus-allele to snp info
    snpinfo = {}
    locusfiles = {}
    loc_to_snppos = {}
    loc_snp_ref = {}

    # build locus allele sql
    orlist = []

    for loc in loc2allls:
        for allele in loc2allls[loc]:
            if allele not in ("1", "0"):
                match = """("{app}_allele".locus_id = '{loc}' AND "{app}_allele".identifier = '{allele}')""".format(
                    app=args.appname, loc=loc, allele=allele)
                orlist.append(match)

    # for ap in ap2all:
    #     for loc in ap2all[ap]:
    #         allele = ap2all[ap][loc]
    #         match = """("{app}_allele".locus_id = '{loc}' AND "{app}_allele".identifier = '{allele}')""".format(app=args.appname,loc=loc,allele=allele)
    #         orlist.append(match)
    orstring = " OR ".join(orlist)

    sqlq = """ 
    SELECT "{a}_allele".locus_id,"{a}_allele".identifier,"{a}_snp".original_aa,"{a}_snp".altered_aa,"{a}_snp".position,"{a}_allele".file_location 
    FROM "{a}_allele" JOIN "{a}_allele_snps" ON "{a}_allele".id = "{a}_allele_snps".allele_id
    JOIN "{a}_snp" ON "{a}_allele_snps".snp_id = "{a}_snp".id
    WHERE ( {orls} );
    """.format(a=args.appname, orls=orstring)

    res = sqlquery_to_outls(conn, sqlq)
    c = 0
    for snp in res:
        loc = snp[0]
        allele = snp[1]
        pos = snp[4]
        ref = snp[2]
        mut = snp[3]
        file = snp[5]

        if loc not in loc_to_snppos:
            loc_to_snppos[loc] = [int(pos)]
        else:
            if pos not in loc_to_snppos[loc]:
                loc_to_snppos[loc].append(int(pos))

        if loc not in loc_snp_ref:
            loc_snp_ref[loc] = {int(pos): ref}
        else:
            if pos not in loc_snp_ref[loc]:
                loc_snp_ref[loc][pos] = ref

        if loc not in locusfiles:
            locusfiles[loc] = file
        if loc not in snpinfo:
            snpinfo[loc] = {allele: {pos: [(ref, mut)]}}
        else:
            if allele not in snpinfo[loc]:
                snpinfo[loc][allele] = {pos: [(ref, mut)]}
            else:
                if pos not in snpinfo[loc][allele]:
                    snpinfo[loc][allele][pos] = [(ref, mut)]
                else:
                    if (ref, mut, pos) not in snpinfo[loc][allele][pos]:
                        snpinfo[loc][allele][pos].append((ref, mut))

    # for loc in snpinfo:
    #     print(loc,snpinfo[loc])
    nloc_to_snppos = {}
    for loc in loc_to_snppos:
        nloc_to_snppos[loc] = list(sorted(loc_to_snppos[loc]))

    return snpinfo, locusfiles, nloc_to_snppos, loc_snp_ref


def make_sorted_list(mgt9loci, loci_pos_dict):
    sortedlist = sorted(mgt9loci, key=lambda x: loci_pos_dict[x][0])
    return sortedlist


def make_snpalign(ap2all, loci_pos_dict, snpinfo, mgt9loci, locusfiles, loc_to_snppos, loc_snp_ref):
    """
    for each apid
    for each locus in ap
    if locus in snpinfo
    get snp pos for that locus from loc_to_snppos
    iterate over pos in order
    if ap allele has snp add snp
    if allele is 0 add n
    if allele neg open allele file and check pos
    else add ref
    :param ap2all:
    :param apidToisolate:
    :param snpinfo:
    :return:
    """
    outstring = {}
    outrefpos = []
    locuslis = []
    c = {}
    overall = 0
    # testloc = "STMMW_23411"
    refseq = ""

    sorted_locils = make_sorted_list(mgt9loci, loci_pos_dict)

    for locus in sorted_locils:
        if locus in snpinfo:
            locuslis.append(locus)
            for pos in loc_to_snppos[locus]:
                overall += 1
                refseq += loc_snp_ref[locus][pos]
                outrefpos += [get_ref_pos(locus, pos, loci_pos_dict)]
                if overall % 100 == 0:
                    print(overall, "positions done")
                # if locus == testloc:
                #     print(snpinfo[locus])
                for ap in list(sorted(ap2all.keys())):
                    if ap not in outstring:
                        outstring[ap] = {}
                        c[ap] = 0
                    apallele = ap2all[ap][locus]
                    # if locus == testloc:
                    #     print(ap,apallele)
                    if "-" in apallele:
                        file = locusfiles[locus]
                        loc_alleles = SeqIO.parse(file, "fasta")
                    snp = ""
                    if apallele in snpinfo[locus]:
                        apsnps = snpinfo[locus][apallele]
                        # if locus == testloc:
                        #     print(apsnps)
                        called = False
                        for snppos in apsnps:
                            if snppos == pos:
                                snp = apsnps[snppos][0][1]
                                called = True
                        if not called:
                            snp = loc_snp_ref[locus][pos]
                    else:
                        if apallele == "1":
                            snp = loc_snp_ref[locus][pos]
                        elif apallele == "0":
                            snp = "N"
                        elif "-" in apallele:
                            alleleid = locus + ":" + apallele
                            for i in loc_alleles:
                                if i.id == alleleid:
                                    seq = str(i.seq)
                                    snp = seq[pos]
                        else:
                            snp = loc_snp_ref[locus][pos]

                    if locus not in outstring[ap]:
                        outstring[ap][locus] = snp
                        c[ap] += 1
                    else:
                        outstring[ap][locus] += snp
                        c[ap] += 1

                    # if locus == testloc:
                    #     print(ap,locus,pos,c[ap],overall,snp)
                    #     sl(0.2)

                    # if c[ap] != overall:
                    #     print(ap,locus,pos,c[ap],overall,snp)
                    #     sl(0.2)

    concatout = {}

    for ap in outstring:
        # print(overall,ap,c[ap])
        concatout[ap] = ""
        for locus in locuslis:
            concatout[ap] += outstring[ap][locus]
            # print(ap,locus,outstring[ap][locus])
            # sl(0.1)
    return concatout, refseq, outrefpos


def get_loci_info(args, conn):
    pqsl_query = """SELECT "identifier","start_location","end_location","orientation" from "{}_locus";""".format(
        args.appname)

    loci_list = sqlquery_to_outls(conn, pqsl_query)

    locidict = {x[0]: x[1:] for x in loci_list}

    return locidict


def get_locils(args, conn, maxlev):
    pqsl_query = """SELECT locus_id from "{}_scheme_loci" WHERE "scheme_id" = 'MGT{}';""".format(args.appname, maxlev)

    loci_list = sqlquery_to_outls(conn, pqsl_query)
    loci_list = [x[0] for x in loci_list]
    return loci_list


def main():
    args = parseargs()

    conn = get_conn()
    # DbConString = "dbname='{0}' host='{1}' port='{2}' user='{3}' password='{4}'".format("salmonella",
    #                                                                                     "0.0.0.0",
    #                                                                                     "5433",
    #                                                                                     "postgres",
    #                                                                                     "P8ppacre!")  ## connection info for Db - assign new user that can only do what is needed for script

    # conn = psycopg2.connect(DbConString)
    # conn.autocommit = True

    dash_nodash, nodash_dash = make_dash_nodash(conn, args)
    tablenos = get_table_nos(conn, args)

    # allele_folder = settings.SUBDIR_ALLELES
    maxscheme = get_max_scheme(conn, args)
    mgt9loci = get_locils(args, conn, maxscheme)
    loci_pos_dict = get_loci_info(args, conn)
    apidToisolate, missing = get_apids(args, conn, maxscheme)

    ap2all, locus2allelels = get_apid_to_alleles(args, conn, apidToisolate, maxscheme, tablenos, nodash_dash)

    snpinfo, locusfiles, loc_to_snppos, loc_snp_ref = get_snpinfo(args, conn, locus2allelels)

    ap_snp_align, refseq, outrefpos = make_snpalign(ap2all, loci_pos_dict, snpinfo, mgt9loci, locusfiles, loc_to_snppos,
                                                    loc_snp_ref)

    if args.fasta_output:
        outf = open(args.outfile, "w")

        outf.write(">Reference\n{}\n".format(refseq))

        for ap in apidToisolate:
            for id in apidToisolate[ap]:
                outf.write(">{}\n{}\n".format(id, ap_snp_align[ap]))
        outf.close()
    else:

        outf = open(args.outfile, "w")
        idls = []
        aplist = []
        for ap in list(sorted(apidToisolate.keys())):
            for id in list(sorted(apidToisolate[ap])):
                aplist.append(ap)
                idls.append(id)
        outf.write("Position\tReference\t{}\n".format("\t".join(idls)))

        for i in range(len(outrefpos)):
            refpos = outrefpos[i]

            refposseq = refseq[i]

            outf.write("{}\t{}".format(refpos, refposseq))
            for ap in aplist:
                strainposseq = ap_snp_align[ap][i]
                outf.write("\t" + strainposseq)
            outf.write("\n")
        outf.close()

    # c = 0
    # for i in locus2allelels:
    #     if len(locus2allelels[i]) > 1:
    #         print(i,locus2allelels[i])
    #         # sl(0.5)
    #         c+=1
    # print(c)
    # for i in ap2all:
    #     for j in ap2all[i]:
    #         print(i,j,ap2all[i][j])
    #         sl(0.3)


if __name__ == '__main__':
    main()