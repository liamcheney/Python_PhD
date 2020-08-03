import argparse
from time import sleep as sl
import pandas as pd

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args

def STs_dist_by_cont(df, geography, MGT_level, min_st_size):
    ST_freq = df[MGT_level].value_counts().to_dict()
    MGT_STs = df[MGT_level].unique().tolist()
    continent_list = df[geography].unique().tolist()
    continent_list_f = []
    for x in continent_list:
        if isinstance(x,float):
            pass
        elif ' ' in x:
            fix = x.replace(' ','_')
            continent_list_f.append(fix)
        else:
            continent_list_f.append(x)
    continent_list = continent_list_f

    con_head = '\t'.join(continent_list)
    # print(f"ST {con_head}", sep='\t')

    for ST in MGT_STs:
        ST_sub = df[df[MGT_level] == ST]
        if ST_sub.shape[0] >= min_st_size:
            continent_dist = ST_sub[geography].value_counts().to_dict()

            # check the cont dists for each ST
            total_st = ST_freq[ST]
            save_list = []
            save_list.append(ST)
            for cont in continent_list:
                if cont in continent_dist.keys():
                    size = str(continent_dist[cont])
                    save_list.append(size)

                    #print out a ST and continent if >70 coverage.
                    if float(int(size) / int(total_st) * 100) > 70.0:
                        print(f"ST{str(ST)} {cont} {size} {total_st} {str(round(int(size)/int(total_st)*100,1))}", sep='\t')
                else:
                    save_list.append('0')
            save_list.append(str(total_st))
            out = '\t'.join(str(v) for v in save_list)
            # print(out, sep='\t')

def main():
    args = parseargs()

    ##find the unique MGT3 STs for a continent
    df = pd.read_csv('/Users/liamcheneyy/Desktop/meta_MGT_isolate_data.txt',sep='\t')

    ##find largest STs for a level, then see the continent distribution
    geography = 'Country'
    MGT_level = 'MGT3'
    min_st_size = 10

    STs_dist_by_cont(df, geography, MGT_level, min_st_size)








if __name__ == '__main__':
    main()
