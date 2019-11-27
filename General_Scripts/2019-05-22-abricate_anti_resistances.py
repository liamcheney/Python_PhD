import argparse
from time import sleep as sl
import glob

def get_all_gene_types(args):

    #go over and find the AMR genes defined across the entire dataset
    #returns a non-redundant list

    total_genes_list = []
    strain_list = []
    for filename in glob.iglob(args.input_folder + '/*sv'):
        infile = open(filename,'r').read().splitlines()
        accession = filename.split('/')[-1].strip("_abricate.csv")

        if accession not in strain_list:
            strain_list.append(accession)

        for line in infile[1:]:
            col = line.split(',')

            if col[4] not in total_genes_list and float(col[8]) >= 90.00 and float(col[9]) >= 90.00:
                total_genes_list.append(col[4])

    return total_genes_list, strain_list

def absence_and_presence(AMR_list, strain_list, args):
    result_dict = {}

    for strain in strain_list:
        result_dict[strain] = {}
        ##add the possible genes for strain
        for element in AMR_list:
            result_dict[strain][element] = ""

        ###open file for strain
        filename = args.input_folder + "/" + strain + '_abricate.csv'
        infile = open(filename,'r').read().splitlines()

        for line in infile[1:]:
            col = line.split(',')

            ##find genes and add True, false == leave blank
            if col[4] in AMR_list and float(col[8]) >= 90.00 and float(col[9]) >= 90.00:
                result_dict[strain][col[4]] = 'TRUE'

    return result_dict

def write_out(result_dict, AMR_list, args):

    with open(args.output_folder + '/abricate_absence_presence.csv','w') as out:
        out.write("Accession" + ',')
        for element in AMR_list:
            out.write(element + ',')
        out.write('\n')

        for strain in result_dict.keys():
            out.write(strain + ',')
            for el in AMR_list:
                out.write(result_dict[strain][el] + ',')
            out.write('\n')

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-i", "--input_folder", required=True)
    parser.add_argument("-o", "--output_folder", required=True)

    args = parser.parse_args()

    return args

def main():
    args = parseargs()
    AMR_list, strain_list = get_all_gene_types(args)
    result_dict = absence_and_presence(AMR_list, strain_list, args)
    write_out(result_dict, AMR_list, args)

if __name__ == '__main__':
    main()
