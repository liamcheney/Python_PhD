

from time import sleep as sl
import argparse
import psycopg2
import sys
from os import path

sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))

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

def neg_to_pos(inallele):
    if "-" in inallele:
        pos = inallele.split("_")[0][1:]
        return pos
    else:
        return inallele

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("id1", help="first strain name or ST (requires -l flag) to compare")
    parser.add_argument("id2", help="second strain name or ST (requires -l flag) to compare")
    parser.add_argument("Database", help="second strain name or ST (requires -l flag) to compare")
    parser.add_argument("-s","--strains", help="Set search for strains not STs (equivalent to comparing MGT9)", action='store_true')
    parser.add_argument("-l", "--level", help="MGT level to search (2-9)",default=9)
    parser.add_argument("-n", "--incnegs", help="Include differences only caused by negative alleles",
                        action='store_true')


    args = parser.parse_args()
    return args

def get_st_id_and_st(args,conn):
    """
    get sts of strain ids for id1 and id2
    :return:
    """

    if args.strains:
        s1 = args.id1
        s2 = args.id2

        hgtquery1 = """ SELECT mgt_id FROM "{var1}_isolate" WHERE "identifier" = '{}'; """.format(s1,var1=args.Database)

        hgtquery2 = """ SELECT mgt_id FROM "{}_isolate" WHERE "identifier" = '{}'; """.format(args.Database, s2)

        res1 = sqlquery_to_outls(conn,hgtquery1)[0][0]

        res2 = sqlquery_to_outls(conn, hgtquery2)[0][0]

        stquery1 = """ SELECT "ap{0}_0","ap{0}_0_st","ap{0}_0_dst" FROM "{1}_view_apcc" WHERE "mgt_id" = '{2}'; """.format(
            str(args.level), args.Database, res1)

        stquery2 = """ SELECT "ap{0}_0","ap{0}_0_st","ap{0}_0_dst" FROM "{1}_view_apcc" WHERE "mgt_id" = '{2}'; """.format(
            str(args.level), args.Database, res2)

        stdst1 = sqlquery_to_outls(conn, stquery1)[0]

        stdst2 = sqlquery_to_outls(conn, stquery2)[0]
    else:
        res1 = tuple(args.id1.split("."))
        res2 = tuple(args.id2.split("."))
        print(res1)
        print(res2)
        stquery1 = """ SELECT "ap{0}_0" FROM "{1}_view_apcc" WHERE ("ap{0}_0_st" = '{2}' AND "ap{0}_0_dst" = {3}); """.format(
            str(args.level), args.Database, res1[0],res1[1])


        stquery2 = """ SELECT "ap{0}_0" FROM "{1}_view_apcc" WHERE ("ap{0}_0_st" = '{2}' AND "ap{0}_0_dst" = {3}); """.format(
            str(args.level), args.Database, res2[0], res2[1])


        try:
            stid1 = sqlquery_to_outls(conn, stquery1)[0][0]
        except:
            print('Error: first ST and dST does not exist')
            sys.exit()

        try:
            stid2 = sqlquery_to_outls(conn, stquery2)[0][0]
        except:
            print('Error: secound ST and dST does not exist')
            sys.exit()

        stdst1 = (stid1,res1[0],res1[1])

        stdst2 = (stid2, res2[0], res2[1])

    return stdst1,stdst2

def get_scheme_tables(args,conn):
    """
    get the allele profile tables as a list that correspond to the correct mgt level
    :param args:
    :return:
    """
    ap_table_query = """ SELECT "table_name" FROM "{0}_tables_ap" WHERE "scheme_id" = 'MGT{1}' """.format(args.Database,args.level)

    tables = sqlquery_to_outls(conn,ap_table_query)

    tables = [x[0] for x in tables]

    tables = sorted(tables, key=lambda n: int(n[-1]))

    return tables

def get_ap(args,conn,st,ap_tables):
    """
    get allele profiles as dict for each locus (key = st.dst)
    :param args:
    :return:
    """

    apid = st[0]
    aplis = []
    locusnames = []

    for table in ap_tables:
        # print(table)
        loci_query = """
        SELECT *
        FROM
        information_schema.columns
        WHERE
        table_name = '{}_{}'
        """.format(args.Database,table)
        locils = sqlquery_to_outls(conn,loci_query)
        locils = [x[3] for x in locils]
        if 'dst' in locils:
            loci = locils[3:-3]
        else:
            loci = locils[1:-6]
        # print(loci[0],loci[-1])

        locusnames += loci
        if "_0" in table:
            allele_query = """
                    SELECT *
                    FROM
                    "{0}_{1}"
                    WHERE "id" = '{2}'
                    """.format(args.Database, table,apid)

            allelels = sqlquery_to_outls(conn, allele_query)
            alleles = allelels[0][3:-3]
            # print(alleles)
            aplis += alleles
        else:
            allele_query = """
                                SELECT *
                                FROM
                                "{0}_{1}"
                                WHERE "main_id" = '{2}'
                                """.format(args.Database, table, apid)

            allelels = sqlquery_to_outls(conn, allele_query)
            alleles = allelels[0][1:-2]
            aplis += alleles

    # print(len(locusnames),len(aplis))

    apdict = dict(zip(locusnames, aplis))

    # for i in apdict:
    #     print(i,apdict[i])
    #     sl(0.5)

    return apdict

def main():
    args = parseargs()

    DbConString = "dbname='vcseventh_11' host='0.0.0.0' port='5432' user='postgres' password='XXXX'"  ## connection info for Db - assign new user that can only do what is needed for script

    conn = psycopg2.connect(DbConString)


    st1,st2 = get_st_id_and_st(args,conn)

    ap_table_list = get_scheme_tables(args,conn)

    apdict1 = get_ap(args,conn,st1,ap_table_list)
    apdict2 = get_ap(args, conn, st2, ap_table_list)

    # print("locus",st1[1:],st2[1:])
    for i in apdict1:
        al1 = apdict1[i]
        al2 = apdict2[i]
        a1pos = str(neg_to_pos(al1))
        a2pos = str(neg_to_pos(al2))
        if args.incnegs:
            if al1 != al2:
                print(i,al1,al2)
                sl(0.1)
        else:
            if a1pos != a2pos:
                print(i,al1,al2)
                sl(0.1)

        # else:
        #     print(i, al1, al2,"match")

if __name__ == '__main__':
	main()



