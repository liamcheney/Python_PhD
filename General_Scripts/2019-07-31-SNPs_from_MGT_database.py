import argparse
from time import sleep as sl

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
def setup_connection(args):

    DbConString = "dbname='{}' host='{}' port='{}' user='{}' password='{}'".format(args.database_name,args.host,args.port,args.psqluser,args.password)  ## connection info for Db - assign new user that can only do what is needed for script
    conn = psycopg2.connect(DbConString)

    return DbConString, conn

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args


def main():
    args = parseargs()


if __name__ == '__main__':
    main()
