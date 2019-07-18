import argparse
import glob
from time import sleep as sl


def check_files(args):
    results_dict = {}
    t3SS_presence_list = []
    t6SS_presence_list = []

    for file in glob.iglob(args.input_folder + '/*sv'):
        accession = file.split('/')[-1].split('_')[0]
        infile = open(file).read().splitlines()

        t6ss_proteins_list = ['vipA/mglA','vipB/mglB','vgrG-2','hcp-2','vasA','vasB','vasC','vasD','vasE','vasF','clpB/vasG','vasH','vasI','vasJ','icmF/vasK']

        type_six_count = 0
        for line in infile[1:]:
            col = []
            if file[-3:] == 'tsv':
                col = line.split('\t')
            elif file[-3:] == 'csv':
                col = line.split(',')


            if float(col[8]) > 75:
                if col[4] in t6ss_proteins_list:
                    type_six_count = type_six_count + 1

            if 'T3SS' in line:
                t3SS_presence_list.append(accession)

            if 'T6SS' in line:
                t6SS_presence_list.append(accession)

        results_dict[accession] = type_six_count

    t6SS_presence_list = set(t6SS_presence_list)
    t3SS_presence_list = set(t3SS_presence_list)

    in_both = []

    for i in t6SS_presence_list:
        if i in t3SS_presence_list:
            in_both.append(i)

    print(len(in_both))
    print(len(t3SS_presence_list))

    return results_dict


# def write_out()

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-i", "--input_folder", required=True,
                        help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args


def main():
    args = parseargs()

    return_dict = check_files(args)

    # for key, value in return_dict.items():
    #     print(key, value)


if __name__ == '__main__':
    main()
