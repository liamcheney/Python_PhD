from time import sleep as sl
import psycopg2
import subprocess
import os
import argparse
import webbrowser

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

def create_list_of_tables(conn, args):

    print("Collecting all allele profiles from the MGT database.")

    wanted_mgt_level = args.wanted_mgt_level

    # create a list of tables for each MGT scheme
    table_numbers_search_string = (""" SELECT table_name FROM "{app}_tables_ap" """).format(app=args.app_name)
    table_numbers_return = sqlquery_to_outls(conn, table_numbers_search_string)
    table_numbers_list = [x[0] for x in table_numbers_return]

    final_table_numbers_list = []
    for table in table_numbers_list:
        table_num = table.split('_')
        if wanted_mgt_level in table_num[0]:
            final_table_numbers_list.append(table)

    return final_table_numbers_list, wanted_mgt_level

def get_column_names(conn, table_numbers_list, args):

    column_headers_list = []
    for table in table_numbers_list:
        table_num = table.split('_')

    # for the found table_name extract the column headers

        # handle the first table to remove extra columns
        if '0' in table_num[1]:
            table_header_string = (""" SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS WHERE TABLE_NAME = N'{app}_{table}' """ ).format(app=args.app_name, table=table)
            column_headers = sqlquery_to_outls(conn, table_header_string)
            column_headers_list.append([x[0] for x in column_headers][1:-6])

        # remove column.id from remaining tables
        else:
            table_header_string = (""" SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS WHERE TABLE_NAME = N'{app}_{table}' """ ).format(app=args.app_name, table=table)
            column_headers = sqlquery_to_outls(conn, table_header_string)
            column_headers_list.append([x[0] for x in column_headers][1:-2])

    column_headers_list = [item for x in column_headers_list for item in x]
    column_headers_list[0] = "#ST"
    column_headers_list[1] = "#DST"
    column_headers_list.insert(0, "#Strain")

    return column_headers_list

def get_strain_names(conn, allele_profiles_list, args):

    strain_hgt_string = (""" SELECT T1.identifier,"{app}_ap9_0".st FROM (SELECT "{app}_isolate".identifier, "{app}_{0}".ap9_0_id FROM "{app}_isolate" LEFT JOIN "{app}_{0}" ON "{app}_isolate".{0}_id::varchar = "{app}_{0}".id::varchar) as T1
         LEFT JOIN "{app}_ap9_0" ON T1.ap9_0_id::varchar = "{app}_ap9_0".id::varchar WHERE "{app}_ap9_0".st is not null """).format(args.mgt_database_version, app=args.app_name)
    strain_hgt_result = sqlquery_to_outls(conn, strain_hgt_string)

    #add information to list from SQL tuples
    #if a select amount of strains desired
    if args.infile:
        allele_select_profiles_dict = {}
        strains_in = list([x for x in open(args.infile, 'r').read().splitlines()])
        for key, value in strain_hgt_result:
            if key in strains_in:
                for line in range(1, len(allele_profiles_list), 1):
                    if value == allele_profiles_list[line][0]:
                        allele_select_profiles_dict[key] = allele_profiles_list[line]

        return allele_select_profiles_dict

    else:
        allele_profiles_dict = {}
        for tup_num in range(0, len(strain_hgt_result), 1):
            for line in range(1, len(allele_profiles_list), 1):
                if strain_hgt_result[tup_num][1] == allele_profiles_list[line][0]:
                    allele_profiles_dict[strain_hgt_result[tup_num][0]] = allele_profiles_list[line]

        return allele_profiles_dict

def combine_allele_profiles(conn, table_numbers_list, column_headers_list, wanted_mgt_level, args):

    #list to store all combined profiles in
    total_list = []

    #go over each table for desired scheme
    alleles_dict = {}
    ap_table_count = 0
    for table in table_numbers_list:
        table_num = table.split('_')
        alleles_dict["table_" + table_num[1]] = {}

        # find the number of cc_cols for table
        cc_col_count = find_number_of_cc_columns(table, conn, args)
        remove_col = cc_col_count + 2 #two is for date create and edit cols

        ap_table_search = (""" SELECT * FROM "{0}_ap{1}_{2}" ORDER BY 1 """).format(args.app_name, str(wanted_mgt_level), str(table_num[1]))
        ap_table_return = sqlquery_to_outls(conn, ap_table_search)

        ##handle first table for extra columns
        if '0' in table_num[1]:

            ap_table = [x[:-remove_col] for x in ap_table_return]
        else:
            ap_table = [x[:-remove_col] for x in ap_table_return]

        for element in ap_table:
            alleles_dict["table_" + table_num[1]][element[0]] = element[1:]

    #create a list of the IDS found in every table dict
    id_list = list([x for x in alleles_dict[list(alleles_dict.keys())[0]]])


    # combine the allele profiles from the individual tables
    for id in id_list:
        alleles_list = []
        for key in alleles_dict.keys():
            for el in alleles_dict[key][id]:
                alleles_list.append(el)
        total_list.append(alleles_list)

    total_list = [column_headers_list] + total_list

    return total_list

def find_number_of_cc_columns(table, conn, args):

    #find the number of column with name = "cc_id"
    cc_col_count = 0
    table_num = table.split('_')

    column_headers_str = (""" SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS WHERE TABLE_NAME = N'{app}_{table}' """).format(app=args.app_name, table=table)
    column_headers_result = sqlquery_to_outls(conn, column_headers_str)
    column_headers = [x[0] for x in column_headers_result]

    #count number of cc_id columns
    for cell in column_headers:
        if "cc" in cell and "_id" in cell:
            cc_col_count = cc_col_count + 1

    return cc_col_count

def convert_to_grapetree(allele_profiles_dict):
    print("Converting allele profiles into Grapetree format.")

    #go over each allele profile
    for key, value in allele_profiles_dict.items():
        for cell_num in range(3, len(allele_profiles_dict[key]), 1):

            #if neg and N allele assignmnts from MGT present then fix for grapetree
            if '_' in allele_profiles_dict[key][cell_num]:
                fix = allele_profiles_dict[key][cell_num].split('_')[0].strip('-')
                allele_profiles_dict[key][cell_num] = fix

            elif 'N' in allele_profiles_dict[key][cell_num]:
                allele_profiles_dict[key][cell_num] = '0'


    return allele_profiles_dict

def generate_grapetree(allele_profiles_dict, column_headers_list, args):
    if not args.skip:
        print("Generating phylogeny for " + str(len(allele_profiles_dict.keys())) + " using Grapetree.")

        ##write out temp file
        with open(args.output_folder + '/allele_profiles.tsv','w') as out:
            for cell in column_headers_list:
                out.write(str(cell) + '\t')
            out.write('\n')

            for key, value in allele_profiles_dict.items():
                out.write(str(key) + '\t')
                for cell in allele_profiles_dict[key]:
                    out.write(str(cell) + '\t')
                out.write('\n')



        grapetree_cmd_string = ("{0} ; {1} grapetree -m {2} -p {3} > {4}").format(args.penv, args.grapetree_path, args.gm, (args.output_folder + "/allele_profiles.tsv"), (args.output_folder + "/grapetree.nwk"))
        print("Using Grapetree command: " + grapetree_cmd_string)
        p = subprocess.Popen(grapetree_cmd_string, shell=True, stdout=subprocess.PIPE)
        out, err = p.communicate()

    else:
        print("Skipping phylogeny generation. Writing out temp file.")

        ##write out temp file
        with open(args.output_folder + '/allele_profiles.tsv','w') as out:
            for cell in column_headers_list:
                out.write(str(cell) + '\t')
            out.write('\n')

            for key, value in allele_profiles_dict.items():
                out.write(str(key) + '\t')
                for cell in allele_profiles_dict[key]:
                    out.write(str(cell) + '\t')
                out.write('\n')

        print("Wrote out allele profiles for "  + str(len(allele_profiles_dict.keys())) + " strains.")

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #Database connecting
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

    #Script running
    parser.add_argument("-m","--wanted_mgt_level", default="9", type=str,
                        help="The MGT level to extract allele profiles and compare alleles to generate phylonegy.")
    parser.add_argument("-o", "--output_folder", required=True,
                        help="Output folder to save phylogeny and alleles profiles.")
    parser.add_argument("-f", "--infile",
                        help="A list of strains accessions to extract from MGT database.")
    parser.add_argument("-dv", "--mgt_database_version", required=True, default="hgt",
                        help="MGT database version. To allow changes in table names to not cause SQL errors.")

    #Grapetree running
    parser.add_argument("-gp_path", "--grapetree_path", default="",
                        help="The path to which grapetree has been installed. If in activate $PATH leave blank")
    parser.add_argument("-skip", default=False,
                        help="Skip the grapetree analysis. Only write out alleles profile.")
    parser.add_argument("-penv", required=True,
                        help="""Command to activate python environment to run Grapetree. Enclose in double quotes. Eg "source activate python_3_env".""")
    parser.add_argument("-gm", default="RapidNJ",
                        help="""Method for grapetre phyloengy model" (Eg. RapidNJ).""")

    args = parser.parse_args()

    return args

def main():

    args = parseargs()

    ##setup the database connection
    DbConString, conn = setup_connection(args)

    #create list of tables
    table_numbers_list, wanted_mgt_level = create_list_of_tables(conn, args)

    #create a list of column headers from all tables
    column_headers_list = get_column_names(conn, table_numbers_list, args)

    #create a list of all allele profiles. multiple tables are combined
    allele_profiles_list = combine_allele_profiles(conn, table_numbers_list, column_headers_list, wanted_mgt_level, args)

    #get strain names
    allele_profiles_dict = get_strain_names(conn, allele_profiles_list, args)

    #convert into grapetree compatible
    allele_profiles_dict = convert_to_grapetree(allele_profiles_dict)

    #create and writeout a nwk phylogeny, allele profile
    generate_grapetree(allele_profiles_dict, column_headers_list, args)

if __name__ == '__main__':
    main()

