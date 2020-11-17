import argparse
from time import sleep as sl


def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args


##will go over a large SRA metadata file. Create a dict of keys per strain. Then output for each isolate.
##essentially converting XML into a TSV.
##no editing or checking for fields

def main():
    args = parseargs()

    infile = open('/Users/liamcheneyy/Downloads/biosample_result.txt').read()

    files = infile.split('\n\n')

    save_dict = {}
    col_heads = []
    for file in files:
        lines = file.split('\n')

        Biosample = lines[1].split(':')[2].split(';')[0].strip(' ')
        save_dict[Biosample] = {}

        #get biosample
        for cell in lines:
            ##handle top
            if '    /' not in cell:
                # print(cell)
                # sl(1)
                pass

           ## handle bottom
            else:
                head = cell.strip('    /').split('=')[0]
                cont = cell.strip('    /').split('=')[1].strip('"')
                save_dict[Biosample][head] = cont
                if head not in col_heads:
                    col_heads.append(head)

    out = open('/Users/liamcheneyy/Downloads/all_SRA_meta.csv','w')

    out.write('BioSample' + ',')
    for el in col_heads:
        out.write(el + ',')
    out.write('\n')

    for key, value in save_dict.items():
        out.write(key + ',')
        for i in col_heads:
            if i in value.keys():
                out.write(save_dict[key][i] + ',')
            else:
                out.write(',')
        out.write('\n')


if __name__ == '__main__':
    main()
