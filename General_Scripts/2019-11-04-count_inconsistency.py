from time import sleep as sl
import argparse
import psycopg2
import matplotlib.pyplot as plt

from statistics import mean
import sys

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    ##Database connecting
    parser.add_argument("-d", "--database_name", required=True,
                        help="sql database to search (eg. vcseventh)")
    parser.add_argument("-a", "--app_name", required=True,
                        help="MGT database App Name(eg. Vibrio)")
    parser.add_argument("-s", "--host",
                        help="db host address", default="0.0.0.0")
    parser.add_argument("-p", "--port",
                        help="psql port", default="5432")
    parser.add_argument("-u", "--psqluser",
                        help="pqsl username", default="postgres")
    parser.add_argument("-w", "--password",
                        help="psql password", required=True)
    parser.add_argument("-n", "--scheme_number",
                        help="Number of scheme levels in MGT database.", default=8)


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

def get_isolate_st_dict(conn, args):

    sqlquery = """ SELECT  "identifier","ap2_0_st","ap3_0_st","ap4_0_st","ap5_0_st","ap6_0_st","ap7_0_st","ap8_0_st","ap{max}_0_st" FROM "{app}_view_apcc" INNER JOIN "{app}_isolate" ON "{app}_view_apcc".mgt_id = "{app}_isolate".mgt_id; """.format(app=args.app_name, max=args.scheme_number)

    res = sqlquery_to_outls(conn,sqlquery)
    isolate_to_mgtst = {}

    for isolatedata in res:
        c = 1
        isolate = isolatedata[0]
        isolate_to_mgtst[isolate] = {}
        for column in isolatedata:
            if c !=1:
                isolate_to_mgtst[isolate][c] = column
            c+=1
    return isolate_to_mgtst

def count_subset_for_lev(lev,isolate_to_mgtst):

    st_to_isolate = {}
    stsize = {}

    matching = {}
    missmatch = {}
    singlemissmatch = {}
    singletons = {}
    zeros = {}
    missmatch_isolates = {}

    if int(lev-1) > 1:
        for prev in range(2,lev):
            matching[prev] = 0
            missmatch[prev] = 0
            singletons[prev] = 0
            zeros[prev] = 0
            singlemissmatch[prev] = 0
            missmatch_isolates[prev] = []
            for isolate in isolate_to_mgtst:
                # print(isolate)
                st = isolate_to_mgtst[isolate][lev]
                if st not in st_to_isolate:
                    st_to_isolate[st] = [isolate]
                else:
                    if isolate not in st_to_isolate[st]:
                        st_to_isolate[st].append(isolate)

            ## currently counting lots of STs twice because some STs that are counted as a minor
            for st in st_to_isolate:
                if st != "0":
                    isolatels = st_to_isolate[st]
                    # print(st,isolatels)
                    # sl(0.5)

                    prev_lev_stlist = {}

                    if len(isolatels) == 1:
                        singletons[prev]+=1
                    else:
                        prev_counted = []
                        for isolate in isolatels:
                            # print(isolate)
                            prev_st = isolate_to_mgtst[isolate][prev]
                            # print(st,prev_st)
                            # sl(2)
                            if prev_st != "0":
                                if prev_st not in prev_lev_stlist:
                                    prev_lev_stlist[prev_st] = 1
                                else:
                                    prev_lev_stlist[prev_st] += 1

                        if len(prev_lev_stlist.keys()) == 0:
                            zeros[prev] +=1
                        else:
                            prevstls = list(prev_lev_stlist.keys())


                            largestkey = max(prev_lev_stlist.keys(), key=(lambda key: prev_lev_stlist[key]))

                            # print(lev,st,prev, prevstls)

                            if len(prevstls) > 1:

                                for prevst in prevstls:
                                    if prevst == largestkey:
                                        matching[prev] += prev_lev_stlist[prevst]
                                        # print("match",prevst,prev_lev_stlist[prevst])
                                    else:
                                        if prev_lev_stlist[prevst] == 1:
                                            missmatch[prev] += 1
                                        else:
                                            missmatch[prev] += prev_lev_stlist[prevst]
                                        for isolate in isolatels:
                                            if prevst == isolate_to_mgtst[isolate][prev]:
                                                missmatch_isolates[prev].append(isolate)


                                        # print("miss",prevst,prev_lev_stlist[prevst])
                                        # print(lev,largestkey,prev_lev_stlist[largestkey],prevst,prev_lev_stlist[prevst])
                                        # sl(2)
                            else:
                                matching[prev] += prev_lev_stlist[largestkey]

                            # sl(0.6)

    return matching,missmatch,singletons,zeros,missmatch_isolates,singlemissmatch

def count_st_size(lev,isolate_to_mgtst,st_sizes):
    st_to_isolate = {}
    st_sizes[lev] = {}

    for isolate in isolate_to_mgtst:
        # print(isolate)
        st = isolate_to_mgtst[isolate][lev]
        if st not in st_to_isolate:
            st_to_isolate[st] = [isolate]
            st_sizes[lev][st] = 1
        else:
            if isolate not in st_to_isolate[st]:
                st_to_isolate[st].append(isolate)
                st_sizes[lev][st] += 1
    return st_sizes

def main():
    """
    for each


    :return:
    """

    args = parseargs()

    DbConString = "dbname='{0}' host='{1}' port='{2}' user='{3}' password='{4}'".format(args.database_name,args.host,args.port,args.psqluser,args.password)  ## connection info for Db - assign new user that can only do what is needed for script

    conn = psycopg2.connect(DbConString)
    conn.autocommit = True

    isolate_to_mgtst = get_isolate_st_dict(conn, args)

    isolate_mismatch_no = {}
    stsizes = {}
    # for lev in range(2, 10):
    #     outf = open("{}_stsizes.txt".format(lev),"w")
    #     stsizes = count_st_size(lev,isolate_to_mgtst,stsizes)
    #     stsizels = list([stsizes[lev][x] for x in stsizes[lev] if stsizes[lev][x]])
    #
    #     print(mean(stsizels),stsizels)
    #     outf.write("\n".join(list(map(str,stsizels))))
    #     outf.close()

    max_level = int(args.scheme_number) - 1

    for lev in range(3,max_level):
        matching,missmatch,singletons,zeros,missmatch_isolates,singlemissmatch = count_subset_for_lev(lev,isolate_to_mgtst)
        sizes = []

        for compare in matching:
            if compare == lev-1:
                total = matching[compare]+missmatch[compare]+singletons[compare]+singlemissmatch[compare]
                missmatch_perc = ((singlemissmatch[compare]+missmatch[compare])/total)*100
                missmatch[compare] =missmatch[compare] + singlemissmatch[compare]

                print(round(missmatch_perc,1))

                print(lev,compare,matching[compare],missmatch[compare],singlemissmatch[compare],singletons[compare])

                for isolate in missmatch_isolates[compare]:
                    pair = '{}_{}'.format(compare,lev)
                    if isolate not in isolate_mismatch_no:
                        isolate_mismatch_no[isolate] = [pair]
                    else:
                        isolate_mismatch_no[isolate].append(pair)

    missmatch_counts = {0:0}

    for isolate in isolate_to_mgtst:
        if isolate in isolate_mismatch_no:
            # nomissmatch = str(",".join(list(sorted(isolate_mismatch_no[isolate]))))
            nomissmatch = len(isolate_mismatch_no[isolate])
            if nomissmatch not in missmatch_counts:
                missmatch_counts[nomissmatch] = 1
            else:
                missmatch_counts[nomissmatch] += 1
        else:
            missmatch_counts[0] += 1

    # missmatch_counts = {x:missmatch_counts[x] for x in missmatch_counts.keys() if len(x.split(","))==1}

    plt.bar(range(len(missmatch_counts)), list(missmatch_counts.values()), align='center')
    plt.xticks(range(len(missmatch_counts)), list(missmatch_counts.keys()))
    plt.savefig("testplot_hier_incocnsistency.pdf")
    print(missmatch_counts)

if __name__ == '__main__':
    main()