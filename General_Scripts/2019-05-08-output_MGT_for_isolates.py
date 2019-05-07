from time import sleep as sl
import argparse
import psycopg2
from Bio import SeqIO
import sys
import itertools
import operator
from os import path
import argparse


def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("database",
                        help="sql database to search (eg salmonella or vcseventh)")
    parser.add_argument("appname",
                        help="application name to search (eg Salmonella or Vibrio)")
    parser.add_argument("-o","--outfile",
                        help="path to file for output")
    parser.add_argument("-s","--host",
                        help="db host address",default="0.0.0.0")
    parser.add_argument("-p","--port",
                        help="psql port",default="5432")
    parser.add_argument("-u","--psqluser",
                        help="pqsl username",default="postgres")
    parser.add_argument("-w","--password",
                        help="psql password",required=True)


    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-f',
                       '--file',
                       help="a file with strain names to retreive one per line. Mutually exclusive to -l,--seqtype,--cc,--odc")
    group.add_argument('-l',
                       '--list',
                       help="A comma separated list of name ids. Mutually exclusive to -f,--seqtype,--cc,--odc")
    group.add_argument('--seqtype',
                       help="sequence type and level separated by comma which will return all matches (i.e. 234.1,8). Mutually exclusive to -l,-f,--cc,--odc")
    group.add_argument('--cc',
                       help="clonal complex and level separated by comma which will return all matches (i.e. 220,7). Mutually exclusive to -l,-f,--seqtype,--odc")
    group.add_argument('--odc',
                       help="odc and odc cutoff separated by comma which will return all matches (i.e. 220,2 or 256,10). Mutually exclusive to -l,-f,--seqtype.--cc")


    args = parser.parse_args()
    return args

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


def mgt_from_ids(args,ids,conn):

    mgtdict = {}

    # print(ids)


    sqlq = """ SELECT "identifier","hgt_id" FROM "{0}_isolate" WHERE "identifier" in ('{1}') """.format(args.appname,"','".join(ids))
    res = sqlquery_to_outls(conn,sqlq)
    for out in res:
        ident = out[0]
        mgt = out[1]
        mgtdict[ident] = mgt

    return mgtdict



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


def mgt_from_assignments(args,conn):

    mgts = []

    if args.seqtype:
        stdst = args.seqtype.split(",")[0]
        level = args.seqtype.split(",")[1]
        st = stdst[0]
        if len(stdst) == 1:
            sqlq = """ SELECT "hgt_id" FROM "{}_view_apcc" WHERE "ap{}_0_st" = '{}'; """.format(args.appname,level,st)
        else:
            dst = stdst[1]
            sqlq = """ SELECT "hgt_id" FROM "{}_view_apcc" WHERE ("ap{}_0_st" = '{}' AND "ap{}_0_st" = '{}')""".format(args.appname,level,st,level,dst)

        res = sqlquery_to_outls(conn,sqlq)

        mgts = [str(x[0]) for x in res]

        ## get st,cc,odcs from "Salmonella_view_apcc" and get strain names from "Salmonella_isolate" by matching mgts

    elif args.cc:
        cc = [str(args.cc.split(",")[0])]
        # print(cc)
        level = str(args.cc.split(",")[1])

        cclis = get_merge_cclis(conn,cc,level,args).split(",")

        # print(cclis)

        sqlq = """ SELECT "hgt_id" FROM "{}_view_apcc" WHERE "cc1_{}" IN ('{}'); """.format(args.appname, level, "','".join(cclis))

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

        sqlq = """ SELECT "hgt_id" FROM "{}_view_apcc" WHERE "cc2_{}" IN ('{}'); """.format(args.appname, odclev, "','".join(odclis))

        res = sqlquery_to_outls(conn, sqlq)

        mgts = [str(x[0]) for x in res]
        ## get st,cc,odcs from "Salmonella_view_apcc" and get strain names from "Salmonella_isolate" by matching mgts


    sqlq = """ SELECT "identifier","hgt_id" FROM "{0}_isolate" WHERE "hgt_id" IN ('{1}') """.format(args.appname, "','".join(mgts))
    out = sqlquery_to_outls(conn, sqlq)

    mgtdict = {}
    ids = []
    for x in out:
        id = x[0]
        mgtid = x[1]
        mgtdict[id] = mgtid
        ids.append(id)

    return ids,mgtdict

def get_mgt_idlist(args,conn):

    if args.file or args.list:
        if args.file:
            ids = []
            with open(args.file,'r') as inf:
                for line in inf:
                    ids.append(line.strip("\n"))
        elif args.list:
            ids = args.list.split(",")
        mgtdict = mgt_from_ids(args,ids,conn)


    elif args.seqtype or args.cc or args.odc:
        ids, mgtdict = mgt_from_assignments(args,conn)

    return ids,mgtdict

# def get_lowest_mergecc(conn,args,cc,level):
#
#     if cc != "None":
#         cc = [cc]
#         ccls = get_merge_cclis(conn,cc,level,args.appname)
#         ccls = ccls.split(",")
#         ccls = [x for x in ccls if x != '0']
#         # print(ccls)
#         ccmin = min(map(int,ccls))
#
#     else:
#         ccmin = "None"
#
#     return ccmin

# def get_lowest_mergeodc(conn,args,cc,level):
#
#     if cc != "None":
#         cc = [cc]
#         ccls = get_merge_odclis(conn,cc,level,args.appname)
#         ccls = ccls.split(",")
#         ccls = [x for x in ccls if x != 'None']
#         # print(ccls)
#         ccmin = min(map(int,ccls))
#
#     else:
#         ccmin = "None"
#
#     return ccmin


def recursive_mergels(ccls,mergels,newmerge):
    if newmerge == 0:
        return ccls
    else:
        nccls = list(ccls)
        new = 0
        for cc in ccls:
            for pair in mergels:
                if cc in pair:
                    if pair[0] == cc:
                        if pair[1] not in nccls:
                            # print(pair)
                            nccls.append(pair[1])
                            new += 1
                    elif pair[1] == cc:
                        if pair[0] not in nccls:
                            # print(pair)
                            nccls.append(pair[0])
                            new += 1
        # print(len(nccls),new)
        return recursive_mergels(nccls,mergels,new)


def get_merges(args,conn):
    # cc merges

    sqlquery = """ SELECT DISTINCT "scheme_id" 
            FROM "{0}_tables_ap";""".format(args.appname)

    res = sqlquery_to_outls(conn, sqlquery)

    schemels = [x[0] for x in res]

    ccmerges = {}

    for scheme in schemels:
        lev = scheme[3:]
        # print("\n\n{}\n\n".format(lev))

        ccmerges[lev] = {}

        sqlquery = """ SELECT DISTINCT "identifier","merge_id_id" 
            FROM "{0}_cc1_{1}";""".format(args.appname,lev)

        res = sqlquery_to_outls(conn, sqlquery)

        mergels = [[x[0],x[1]] for x in res if len(x) > 1 and x[1]]
        # print(mergels)

        for ccpair in mergels:
            cc1 = ccpair[0]
            cc2 = ccpair[1]
            if cc1 not in ccmerges[lev] and cc2 not in ccmerges[lev]:
                ccAllMerges = recursive_mergels([cc1],mergels,1)
                for i in ccAllMerges:
                    i=str(i)
                    ccmerges[lev][i] = ccAllMerges


    # odc merges

    odcmerges = {}

    sqlquery = """ SELECT DISTINCT "table_name" 
            FROM "{0}_tables_cc" WHERE SUBSTRING("table_name",3,1) = '2';""".format(args.appname)

    res = sqlquery_to_outls(conn, sqlquery)

    tables = [x[0] for x in res]

    for tab in tables:
        # print("\n\n{}\n\n".format(tab))

        odcmerges[tab] = {}
        sqlquery = """ SELECT DISTINCT "identifier","merge_id_id" 
            FROM "{0}_{1}";""".format(args.appname,tab)

        res = sqlquery_to_outls(conn, sqlquery)

        mergels = [[x[0],x[1]] for x in res if len(x) > 1 and x[1]]

        for odcpair in mergels:
            odc1 = odcpair[0]
            odc2 = odcpair[1]
            if odc1 not in odcmerges[tab] and odc2 not in odcmerges[tab]:
                odcAllMerges = recursive_mergels([odc1],mergels,1)
                for i in odcAllMerges:
                    i = str(i)
                    odcmerges[tab][i] = odcAllMerges

    return ccmerges,odcmerges



def get_info_writeout(args,conn,ids,mgtdict,ccmerges,odcmerges):
    loci_query = """
    SELECT *
    FROM
    information_schema.columns
    WHERE
    table_name = '{}_view_apcc'
    """.format(args.appname)

    res = sqlquery_to_outls(conn,loci_query)

    sts = [x[3] for x in res if "_st" in x[3]]
    dsts = [x[3] for x in res if "_dst" in x[3]]
    ccs = [x[3] for x in res if "cc1_" in x[3] and "merge" not in x[3]]
    odcs = [x[3] for x in res if "cc2_" in x[3]and "merge" not in x[3]]

    # print(sts)
    # print(dsts)
    # print(ccs)
    # print(odcs)

    mgt_to_id = {}
    for k, v in mgtdict.items():
        mgt_to_id.setdefault(v, []).append(k)

    mgtls = mgt_to_id.keys()

    sts_sqlget = """ SELECT "hgt_id","{0}" FROM "{1}_view_apcc" WHERE "hgt_id" IN ('{2}') """.format('","'.join(map(str,sts)),args.appname,"','".join(map(str,mgtls)))

    stres = sqlquery_to_outls(conn,sts_sqlget)

    dsts_sqlget = """ SELECT "hgt_id","{0}" FROM "{1}_view_apcc" WHERE "hgt_id" IN ('{2}') """.format('","'.join(map(str,dsts)),args.appname,"','".join(map(str,mgtls)))

    dstres = sqlquery_to_outls(conn,dsts_sqlget)

    ccs_sqlget = """ SELECT "hgt_id","{0}" FROM "{1}_view_apcc" WHERE "hgt_id" IN ('{2}') """.format('","'.join(map(str,ccs)),args.appname,"','".join(map(str,mgtls)))

    ccres = sqlquery_to_outls(conn,ccs_sqlget)

    odcs_sqlget = """ SELECT "hgt_id","{0}" FROM "{1}_view_apcc" WHERE "hgt_id" IN ('{2}') """.format('","'.join(map(str,odcs)),args.appname,"','".join(map(str,mgtls)))

    odccres = sqlquery_to_outls(conn,odcs_sqlget)

    mgt_sts = {x[0]:x[1:] for x in stres}

    mgt_ccs = {x[0]:x[1:] for x in ccres}

    mgt_odcs = {x[0]:x[1:] for x in odccres}

    mgt_stdsts = {}

    for i in zip(stres,dstres):
        mgt = ''
        c=0
        for j in zip(i[0],i[1]):
            if c == 0:
                mgt = j[0]
                mgt_stdsts[mgt] = []
                c+=1
            else:
                stdst = str(j[0])+"."+str(j[1])
                mgt_stdsts[mgt].append(stdst)

    # TODO get smallest merged for each level of cc output

    stheads = ["MGT2 ST", "MGT3 ST", "MGT4 ST", "MGT5 ST", "MGT6 ST", "MGT7 ST", "MGT8 ST", "MGT9 ST"]
    stdstheads = [x.replace("ST","STDST") for x in stheads]
    ccheads = [x.replace("ST","CC") for x in stheads]
    ccheadsmin = [x.replace("ST", "CC mergemin") for x in stheads]
    odcheads = ["ODC1","ODC2","ODC5","ODC10"]
    odcheadsmin = ["ODC1 mergemin", "ODC2 mergemin", "ODC5 mergemin", "ODC10 mergemin"]

    if args.outfile:
        outf = open(args.outfile,"w")
        outf.write("Strain ID\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format("\t".join(stheads),"\t".join(stdstheads),"\t".join(ccheads),"\t".join(ccheadsmin),"\t".join(odcheads),"\t".join(odcheadsmin)))
        for strain in mgtdict:
            print(strain)
            mgtid = mgtdict[strain]
            stlis = list(map(str,mgt_sts[mgtid]))
            dstlis = list(map(str,mgt_stdsts[mgtid]))
            cclis = list(map(str,mgt_ccs[mgtid]))
            odclis = list(map(str,mgt_odcs[mgtid]))
            ccminls = []
            odcminls = []
            res = {}
            for pos,cc in enumerate(cclis):
                lev = pos+2
                if (lev,cc) in res:
                    ncc = res[(lev,cc)]
                    ccminls.append(str(ncc))
                else:
                    # for i in ccmerges[str(lev)]:
                    #     print(str(lev),i,ccmerges[str(lev)][i])
                    # print(ccmerges[str(lev)])
                    if cc in ccmerges[str(lev)]:
                        ncc = min(map(int,ccmerges[str(lev)][cc]))
                        print(lev,cc,ncc)
                    else:
                        ncc = cc
                    # ncc = get_lowest_mergecc(conn,args,cc,lev)
                    ccminls.append(str(ncc))
                    res[(lev,cc)] = ncc
            res2 = {}



            for pos,odc in enumerate(odclis):
                if pos == 0:
                    continue
                lev = str(pos+1)
                convback = {"1": "1", "2": "2", "3": "5", "4": "10"}
                # lev = convback[lev]
                if (lev,odc) in res2:
                    nodc = res2[(lev,odc)]
                    odcminls.append(str(nodc))
                else:
                    if odc in odcmerges["cc2_"+lev]:
                        nodc = min(map(int, odcmerges["cc2_"+lev][odc]))
                        print(lev,odc,nodc)
                    else:
                        nodc = odc
                    # nodc = get_lowest_mergeodc(conn,args,odc,lev)
                    odcminls.append(str(nodc))
                    res2[(lev,odc)] = nodc




            outf.write(
                "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{7}\t{6}\n".format(strain,"\t".join(stlis), "\t".join(dstlis), "\t".join(cclis), "\t".join(ccminls), "\t".join(odclis),"\t".join(odcminls),ccminls[-1]))
        outf.close()
    else:
        sys.stdout.write("Strain ID\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format("\t".join(stheads),"\t".join(stdstheads),"\t".join(ccheads),"\t".join(ccheadsmin),"\t".join(odcheads),"\t".join(odcheadsmin)))


        for strain in mgtdict:
            print(strain)
            mgtid = mgtdict[strain]
            stlis = list(map(str,mgt_sts[mgtid]))
            dstlis = list(map(str,mgt_stdsts[mgtid]))
            cclis = list(map(str,mgt_ccs[mgtid]))
            odclis = list(map(str,mgt_odcs[mgtid]))
            ccminls = []
            odcminls = []
            res = {}
            for pos, cc in enumerate(cclis):
                lev = pos + 2
                if (lev, cc) in res:
                    ncc = res[(lev, cc)]
                    ccminls.append(str(ncc))
                else:
                    # for i in ccmerges[str(lev)]:
                    #     print(str(lev),i,ccmerges[str(lev)][i])
                    # print(ccmerges[str(lev)])
                    if cc in ccmerges[str(lev)]:
                        ncc = min(map(int, ccmerges[str(lev)][cc]))
                        print(lev, cc, ncc)
                    else:
                        ncc = cc
                    # ncc = get_lowest_mergecc(conn,args,cc,lev)
                    ccminls.append(str(ncc))
                    res[(lev, cc)] = ncc
            res2 = {}

            for pos, odc in enumerate(odclis):
                if pos == 0:
                    continue
                lev = str(pos + 1)
                convback = {"1": "1", "2": "2", "3": "5", "4": "10"}
                # lev = convback[lev]
                if (lev, odc) in res2:
                    nodc = res2[(lev, odc)]
                    odcminls.append(str(nodc))
                else:
                    if odc in odcmerges["cc2_" + lev]:
                        nodc = min(map(int, odcmerges["cc2_" + lev][odc]))
                        print(lev, odc, nodc)
                    else:
                        nodc = odc
                    # nodc = get_lowest_mergeodc(conn,args,odc,lev)
                    odcminls.append(str(nodc))
                    res2[(lev, odc)] = nodc

            sys.stdout.write(
                "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{7}\t{6}\n".format(strain, "\t".join(stlis), "\t".join(dstlis),
                                                             "\t".join(cclis), "\t".join(ccminls), "\t".join(odclis),
                                                             "\t".join(odcminls),ccminls[-1]))



    # print(stres)
    # print(dstres)
    # print(mgt_stdsts)
    # print(ccres)
    # print(odccres)



def main():
    """
    1 - get mgt ids for isolates in input
    2 - use mgt id to get st,dst,cc and odcs
    3 - check merges for lowest merge cc
    return csv
    :return:
    """
    args = parseargs()

    DbConString = "dbname='{}' host='{}' port='{}' user='{}' password='{}'".format(args.database,args.host,args.port,args.psqluser,args.password)  ## connection info for Db - assign new user that can only do what is needed for script

    conn = psycopg2.connect(DbConString)

    ccmerges,odcmerges = get_merges(args,conn)

    ids, mgtdict = get_mgt_idlist(args,conn)

    get_info_writeout(args,conn,ids,mgtdict,ccmerges,odcmerges)

    # print(ids)
    # print(mgtdict)






if __name__ == '__main__':
    main()