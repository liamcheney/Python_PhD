from time import sleep as sl
import argparse
import psycopg2, getpass



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

def get_scheme_tables(args,conn):
    """
    get the allele profile tables as a list that correspond to the correct mgt level
    :param args:
    :return:
    """
    ap_table_query = """ SELECT "table_name" FROM "{0}_tables_ap" WHERE "scheme_id" = 'MGT{1}' """.format(args.appname,args.lev)

    tables = sqlquery_to_outls(conn,ap_table_query)

    tables = [x[0] for x in tables]

    tables = sorted(tables, key=lambda n: int(n[-1]))

    return tables

def get_ap(args,conn,apidlist,ap_tables,maxlevel):
    """
    get allele profiles as dict for each locus (key = st.dst)
    :param args:
    :return:
    """

    aplis = []
    locusnames = []
    aplisdict = {}

    apidlist = list(map(str,apidlist))
    #### TODO the below "X" in table is the largest scheme number
    for table in ap_tables:
        if maxlevel in table:
            last = -6
        else:
            last = -3
        # print(table)
        loci_query = """
        SELECT *
        FROM
        information_schema.columns
        WHERE
        table_name = '{}_{}'
        """.format(args.appname,table)
        locils = sqlquery_to_outls(conn,loci_query)
        locils = [x[3] for x in locils]
        if 'dst' in locils:
            loci = locils[3:last]
        else:
            loci = locils[1:-2]
        # print(loci[0],loci[-1])

        locusnames += loci
        if "_0" in table:
            allele_query = """
                    SELECT *
                    FROM
                    "{0}_{1}"
                    WHERE "id" in ('{2}');
                    """.format(args.appname, table,"','".join(apidlist))

            allelels = sqlquery_to_outls(conn, allele_query)
            for res in allelels:
                res_ap = res[3:last]
                resid = res[0]
                aplisdict[resid] = list(res_ap)

                # alleles = allelels[0][3:last]
            # print(alleles)
            # aplis += alleles
        else:
            allele_query = """
                                SELECT *
                                FROM
                                "{0}_{1}"
                                WHERE "main_id" in ('{2}');
                                """.format(args.appname, table, "','".join(apidlist))

            allelels = sqlquery_to_outls(conn, allele_query)
            for res in allelels:
                res_ap = list(res[1:-2])
                resid = res[0]
                aplisdict[resid]+=res_ap

    # print(len(locusnames),len(aplis))
    apdict = {}
    for id in aplisdict:

        apdict[id] = dict(zip(locusnames, aplisdict[id]))

    # apdict = dict(zip(locusnames, aplis))

    # for i in apdict:
    #     print(i,apdict[i])
    #     sl(0.5)

    return apdict,locusnames

def getSts(strainlist,args,conn):
    """

    :param slist:
    :param args:
    :param conn:
    :return:
    """
    if args.all:
        hgtquery1 = """ SELECT "{appname}_isolate".identifier,"ap{lev}_0","ap{lev}_0_st","ap{lev}_0_dst" FROM "{appname}_isolate" JOIN "{appname}_view_apcc" ON "{appname}_isolate".mgt_id = "{appname}_view_apcc".mgt_id; """.format(
            lev=args.lev, appname=args.appname)
    else:
        hgtquery1 = """ SELECT "{appname}_isolate".identifier,"ap{lev}_0","ap{lev}_0_st","ap{lev}_0_dst" FROM "{appname}_isolate" JOIN "{appname}_view_apcc" ON "{appname}_isolate".mgt_id = "{appname}_view_apcc".mgt_id WHERE "{appname}_isolate".identifier in ('{idlist}'); """.format(lev=args.lev,appname=args.appname, idlist="','".join(strainlist))

    res = sqlquery_to_outls(conn, hgtquery1)

    return res

def countZeroAndNeg(apls):
    z = 0
    n = 0
    for i in apls:
        if i==0:
            z+=1
        elif "-" in i:
            n+=1
    return z,n

def get_max_scheme(connection, args):
    """
    :param connection: sql connection
    :param args: inputs
    :return: highest MGT level in current database
    """
    sqlquery = """ SELECT MAX("display_order") FROM "{}_tables_ap";  """.format(args.appname)

    res = sqlquery_to_outls(connection, sqlquery)

    maxno = int(res[0][0])

    return str(maxno)


def convPos(allele):
    if "-" in allele:
        return allele[1:].split("_")[0]
    else:
        return allele

def get_merge_odclis(connection, odc, odclev,db):
    # print(cc)

    conv = {"1": '1', "2": '2', "5": '3', "10": '4'}

    # odclev = conv[odclev]

    sqlquery = """ SELECT DISTINCT "cc2_{1}","cc2_{1}_merge" 
    FROM "{0}_view_apcc"
    WHERE ("cc2_{1}" in ('{2}') OR "cc2_{1}_merge" in ('{2}'));""".format(db, odclev, "','".join(odc))

    tuplelis = sqlquery_to_outls(connection, sqlquery)
    # print("tuplelis", tuplelis)
    tmpls = [x[0] for x in tuplelis] + [x[1] for x in tuplelis]
    # print("tmppls1",tmpls)
    tmpls = list(set(map(str, tmpls)))
    # print("tmppls2", tmpls)

    tmpls = [x for x in tmpls if x != "None"]

    # print("tmppls3", tmpls)
    cclis = []
    for i in tmpls:
        if i not in cclis:
            cclis.append(i)

    # print(cclis)

    sqlquery = """ SELECT DISTINCT "cc2_{1}","cc2_{1}_merge" 
    FROM "{0}_view_apcc" 
    WHERE ("cc2_{1}" in ('{2}') OR "cc2_{1}_merge" in ('{2}'));""".format(db, odclev, "','".join(cclis))

    tuplelis = sqlquery_to_outls(connection, sqlquery)
    # print("tuplelis", tuplelis)
    tmpls = [x[0] for x in tuplelis] + [x[1] for x in tuplelis]
    # print("tmppls1",tmpls)
    tmpls = list(set(map(str, tmpls)))
    # print("tmppls2", tmpls)

    tmpls = [x for x in tmpls if x != "None"]

    # print("tmppls3", tmpls)
    cclis = []
    for i in tmpls:
        if i not in cclis:
            cclis.append(i)

    # print(cclis)

    cclis = ",".join(cclis)

    return cclis

def id_from_assignments(args,conn):

    mgts = []

    if args.seqtype:
        stdst = args.seqtype.split(",")[0]
        level = args.seqtype.split(",")[1]
        st = stdst[0]
        if len(stdst) == 1:
            sqlq = """ SELECT "mgt_id" FROM "{}_view_apcc" WHERE "ap{}_0_st" = '{}'; """.format(args.appname,level,st)
        else:
            dst = stdst[1]
            sqlq = """ SELECT "mgt_id" FROM "{}_view_apcc" WHERE ("ap{}_0_st" = '{}' AND "ap{}_0_st" = '{}')""".format(args.appname,level,st,level,dst)

        res = sqlquery_to_outls(conn,sqlq)

        mgts = [str(x[0]) for x in res]

        ## get st,cc,odcs from "Salmonella_view_apcc" and get strain names from "Salmonella_isolate" by matching mgts

    elif args.cc:
        cc = [str(args.cc.split(",")[0])]
        # print(cc)
        level = str(args.cc.split(",")[1])

        cclis = get_merge_cclis(conn,cc,level,args).split(",")

        # print(cclis)

        sqlq = """ SELECT "mgt_id" FROM "{}_view_apcc" WHERE "cc1_{}" IN ('{}'); """.format(args.appname, level, "','".join(cclis))

        res = sqlquery_to_outls(conn, sqlq)

        mgts = [str(x[0]) for x in res]
        ## get st,cc,odcs from "Salmonella_view_apcc" and get strain names from "Salmonella_isolate" by matching mgts

    elif args.odc:
        odc = [args.odc.split(",")[0]]
        odclev = args.odc.split(",")[1]
        conv = {"1": '1', "2": '2', "5": '3', "10": '4'}

        odclev = conv[odclev]

        odclis = get_merge_odclis(conn,odc, odclev,args.appname).split(",")

        # print(odclis)

        sqlq = """ SELECT "mgt_id" FROM "{}_view_apcc" WHERE "cc2_{}" IN ('{}'); """.format(args.appname, odclev, "','".join(odclis))

        res = sqlquery_to_outls(conn, sqlq)

        mgts = [str(x[0]) for x in res]
        ## get st,cc,odcs from "Salmonella_view_apcc" and get strain names from "Salmonella_isolate" by matching mgts


    sqlq = """ SELECT "id","identifier","file_alleles" FROM "{0}_isolate" WHERE "mgt_id" IN ('{1}') """.format(args.appname, "','".join(mgts))
    out = sqlquery_to_outls(conn, sqlq)

    straindict = {}
    ids = []
    for x in out:
        id = x[0]
        strainname = x[1]
        alleles = x[2]
        straindict[strainname] = [id,alleles]
        ids.append(strainname)

    return ids,straindict

def get_isolate_idlist(args,conn):

    if args.file or args.list:
        if args.file:
            ids = []
            with open(args.file,'r') as inf:
                for line in inf:
                    ids.append(line.strip("\n"))
        elif args.list:
            ids = args.list.split(",")

        sqlq = """ SELECT "id","identifier","file_alleles" FROM "{0}_isolate" WHERE "identifier" IN ('{1}') """.format(
            args.appname, "','".join(ids))

        out = sqlquery_to_outls(conn, sqlq)

        straindict = {}
        ids = []
        for x in out:
            id = x[0]
            strainname = x[1]
            alleles = x[2]
            straindict[strainname] = [id, alleles]
            ids.append(strainname)

    elif args.seqtype or args.cc or args.odc:
        ids, straindict = id_from_assignments(args,conn)

    elif args.project:
        ids, straindict = id_from_prj(args,conn)

    return ids,straindict

def id_from_prj(args,conn):

    # print(ids)


    sqlq = """ SELECT "{0}_isolate".id,"{0}_isolate".identifier,"{0}_isolate".file_alleles FROM "{0}_isolate" INNER JOIN "{0}_project" ON "{0}_project".id="{0}_isolate".project_id WHERE "{0}_project".identifier = '{1}';""".format(args.appname,args.project)
    out = sqlquery_to_outls(conn,sqlq)

    straindict = {}
    ids = []
    for x in out:
        id = x[0]
        strainname = x[1]
        alleles = x[2]
        straindict[strainname] = [id,alleles]
        ids.append(strainname)

    return ids,straindict

def get_merge_cclis(connection, cc, level,args):
    """
    Gather all CCs that merge with the input CC then gather all CCs that merge with those CCs and return as list
    :param connection: sql connection
    :param cc: input clonal complex to check for merges
    :param level: MGT level
    :param db: database name (= args.psql)
    :return: comma joined list of clonal complexes of original and merges of input cc
    """
    #  Round one cc retreival - get any cc that input cc merges to or that merge to the input cc and store in list
    sqlquery = """ SELECT DISTINCT "identifier","merge_id_id" FROM "{0}_cc1_{1}" WHERE ("identifier" in ('{2}') OR "merge_id_id" in ('{2}'));""".format(args.appname, level, "','".join(cc))

    tuplelis = sqlquery_to_outls(connection, sqlquery)
    # print(cc)
    # print(tuplelis)

    tmpls = [x[0] for x in tuplelis] + [x[1] for x in tuplelis]
    tmpls = list(set(map(str, tmpls)))
    tmpls = [x for x in tmpls if x != "None"]
    tmpls = [x for x in tmpls if x != ""]



    cclis = []
    for i in tmpls:
        if i not in cclis:
            cclis.append(i)
    # print(cclis)
    #  Round two cc retreival - get any cc that round 1 ccs merges to or that merge to the round 1 ccs and store in list
    sqlquery = """ SELECT DISTINCT "identifier","merge_id_id" 
    FROM "{0}_cc1_{1}" 
    WHERE ("identifier" in ('{2}') OR "merge_id_id" in ('{2}'));""".format(args.appname, level, "','".join(cclis))

    tuplelis = sqlquery_to_outls(connection, sqlquery)
    tmpls = [x[0] for x in tuplelis] + [x[1] for x in tuplelis]
    tmpls = list(set(map(str, tmpls)))
    tmpls = [x for x in tmpls if x != "None"]
    tmpls = [x for x in tmpls if x != ""]

    cclis = []
    for i in tmpls:
        if i not in cclis:
            cclis.append(i)
    # print(cclis)
    cclis = ",".join(cclis)

    return cclis

def get_sts_from_stids(args,conn,apids):
    sqlq = """ SELECT "ap{lev}_0","ap{lev}_0_st","ap{lev}_0_dst" FROM "{appname}_view_apcc"  WHERE "ap{lev}_0" in ('{stl}');""".format(lev=args.lev,appname=args.appname, stl="','".join(map(str,apids)))
    res = sqlquery_to_outls(conn, sqlq)
    return res

def get_stids_from_sts(args,conn,sts):

    sqlq = """ SELECT "ap{lev}_0","ap{lev}_0_st","ap{lev}_0_dst" FROM "{appname}_view_apcc"  WHERE "ap{lev}_0_st" in ('{stl}');""".format(lev=args.lev,appname=args.appname, stl="','".join(map(str,sts)))

    res = sqlquery_to_outls(conn, sqlq)
    intact_st = {}
    non_intact_st = {}
    for result in res:
        if result[2] == "0":
            intact_st[result[1]] = result[0]
    for result in res:
        st = result[1]
        stid = result[0]
        if st not in intact_st:
            if st not in non_intact_st:
                non_intact_st[st] = [stid]
            else:
                non_intact_st[st].append(stid)

    return intact_st,non_intact_st


def consolidate_prof_dicts(profdicts):

    outd = {}
    for apid in profdicts:
        profdict = profdicts[apid]
        for locus in profdict:
            if locus not in outd:
                outd[locus] = profdict[locus]
            else:
                if outd[locus] == '0' and profdict[locus] != '0':
                    # print(locus)
                    # print(outd[locus])
                    outd[locus] = profdict[locus]
                    # print("change")
                    # print(profdict[locus],"\n")
    return outd

def consolidateProfs(args,conn,apidls):
    """
    consolidated allele profile:
    25.1 = 1 3 4 2 0 1
    25.2 = 1 0 4 2 8 1
    consolidated ST 25 = 1 3 4 2 8 1

    get consolidated allele profile for each stid in outdict

    1 - get all st.dsts for st that matches stid
    2 - iterate over allele profiles, if non 0 is in location of current 0 then replace
    3 - match new allele profile to stids that have the correct st (new outdict)
    (4) - change dst in non grapetree out to 'c'

    """
    # print(1)
    reslists = get_sts_from_stids(args,conn,apidls)
    # print(2)
    sts = [res[1] for res in reslists]
    sts = list(set(sts))
    intact_st,non_intact_st = get_stids_from_sts(args,conn,sts)
    # print(3)
    outd = {}
    # print(intact_st)
    for res in reslists:
        stid = res[0]
        st = res[1]
        if st in intact_st:
            apid = [intact_st[st]]
            apdict, locusnames = get_ap(args, conn, apid, args.aptables, args.maxlevel)
            apdict = apdict[intact_st[st]]
            outd[stid] = apdict
    # print(4)
    stidls = []
    for st in non_intact_st:
        ids = non_intact_st[st]
        stidls += ids
    apdict, locusnames = get_ap(args, conn, stidls, args.aptables, args.maxlevel)

    for st in non_intact_st:
        st_speficif_apids = {}
        stidls = non_intact_st[st]
        # print(stidls)
        for stid in stidls:
            st_speficif_apids[stid] = apdict[stid]
        cons_apdict = consolidate_prof_dicts(st_speficif_apids)
        for apid in stidls:
            outd[apid] = cons_apdict
    # print(5)
    return outd,locusnames

def parseargs():

    class Password(argparse.Action):
        def __call__(self, parser, namespace, values, option_string):
            if values is None:
                values = getpass.getpass()

            setattr(namespace, self.dest, values)

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # parser.add_argument("--isolatelist", help="one per line list of isolates")
    # /Users/michaelpayne/Desktop/aus_snp_align_test_isolates.txt
    parser.add_argument("--appname", help="App name", default="Staphylococcus")
    parser.add_argument("--dbname", help="database name ", default="saureus_1")
    parser.add_argument("--lev", help="MGT level to get profile for",default="9")
    parser.add_argument("-o", "--outfile", help="output file")
    parser.add_argument("-r", "--port", help="postgres port", default='5432')
    parser.add_argument("-u", "--user", help="postgres user", default='postgres')
    parser.add_argument('-p', action=Password, nargs='?', dest='password', help='Enter your password', default='5678')
    parser.add_argument("--show_nost", help="show strains with no ST at this level", action='store_true')
    parser.add_argument("--consolidated_sts", help="this will make the most complete allele profile for each ST by combining all degenerate STs (where no positive ST exists)", action='store_true')
    parser.add_argument('-g', '--grapetree', help="output format usable for tree building in grapetree",
                        action='store_true')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-f',
                       '--file',
                       help="a file with strain names to retreive one per line. Mutually exclusive to -l,--seqtype,--cc,--odc,--project")
    group.add_argument('-l',
                       '--list',
                       help="A comma separated list of name ids. Mutually exclusive to -f,--seqtype,--cc,--odc,--project")
    group.add_argument('--seqtype',
                       help="sequence type and level separated by comma which will return all matches (i.e. 234.1,8). Mutually exclusive to -l,-f,--cc,--odc,--project")
    group.add_argument('--cc',
                       help="clonal complex and level separated by comma which will return all matches (i.e. 220,7). Mutually exclusive to -l,-f,--seqtype,--odc,--project")
    group.add_argument('--odc',
                       help="odc and odc cutoff separated by comma which will return all matches (i.e. 220,2 or 256,10). Mutually exclusive to -l,-f,--seqtype.--cc,--project")
    group.add_argument('--project',
                       help="project id to return data. Mutually exclusive to -l,-f,--seqtype.--cc")
    group.add_argument('--all',
                       help="get all data in db. Mutually exclusive to -l,-f,--seqtype.--cc")

    args = parser.parse_args()
    return args



def main():
    args = parseargs()

    DbConString = "dbname={} host='0.0.0.0' port='{}' user='{}' password='{}'".format(args.dbname,args.port,args.user,args.password)  ## connection info for Db - assign new user that can only do what is needed for script

    conn = psycopg2.connect(DbConString)
    print("start")
    aptables = get_scheme_tables(args,conn)
    args.aptables = aptables

    maxlevel = get_max_scheme(conn,args)
    args.maxlevel = maxlevel
    strainlist,straindict = get_isolate_idlist(args,conn)

    # if args.isolatelist:
    #     strainlist = open(args.isolatelist,"r").read().splitlines()
    # else:
    #     strainlist = []

    reslists = getSts(strainlist,args,conn)

    print("got STS")
    outdict = {}
    locils = []
    c = 0
    idlist = []
    for inf in reslists:
        stid = inf[1]
        if stid and stid not in idlist:
            idlist.append(stid)
    nost_strains = []

    if args.consolidated_sts:
        print("consolidating")
        outdict,locils = consolidateProfs(args,conn,idlist)
    else:
        outdict, locils = get_ap(args, conn, idlist, aptables, maxlevel)


    if not args.grapetree:
        outf = open(args.outfile,"w")
        outf.write("Strain\tst\tdst\t{}\n".format("\t".join(locils)))

        for res in reslists:
            stid = res[1]
            if stid != None:
                apls = [outdict[stid][x] for x in locils]
                    # z,n = countZeroAndNeg(apls)
                    # print(res[0],res[2],res[3],z,n)
                outf.write("{}\t{}\t{}\t{}\n".format(res[0],res[2],res[3],"\t".join(apls)))
            else:
                    nost_strains.append(res[0])
        outf.close()
    else:
        outf = open(args.outfile,"w")
        outf.write("#Name\t#ST\t{}\n".format("\t".join(locils)))
        for res in reslists:
            stid = res[1]
            if stid != None:
                apls = [outdict[stid][x] for x in locils]
                posapls = [convPos(x) for x in apls]
                # z,n = countZeroAndNeg(apls)
                # print(res[0],res[2],res[3],z,n)
                outf.write("{}\t{}\t{}\n".format(res[0],res[2],"\t".join(posapls)))
            else:
                nost_strains.append(res[0])
        outf.close()
    print("{} strains have no MGT{} ST".format(len(nost_strains),args.lev))
    print("\n".join(nost_strains))


if __name__ == '__main__':
    main()