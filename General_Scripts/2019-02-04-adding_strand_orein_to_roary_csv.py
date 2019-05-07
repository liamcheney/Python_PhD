from time import sleep as sl
import json
import csv
##script to take strand orien for the core genes in each roary core genes group, then add to roary csv which does not have strain oreintation.
##1. make dict with structure dict{group_name[genome_id_gene_id][strand_orein])
##note. will skip falsely split paralogs, since those core genes are excluded from further MSA compliling for PAML analysis. Will save time only doing for genes used in PAML DnDs analysis.
##2. iterate over roary v. cholerae core genome csv and add strand oreintation to each core gene for each core gene group

input_old_csv = '/Users/liamcheneyy/Desktop/D2C_true_core_gene_presence_absence.csv'
dict_out_file = '/Users/liamcheneyy/Desktop/dict.json'
input_cg_csv = "/Volumes/Liam's HDD/PhD/1_core_genomes/4All_Core_Genome/D2C_i96_99_edited_names_percen_para_ana.csv"
outfile_path = '/Users/liamcheneyy/Desktop/D2C_i96_99_new_edited_names_percen_para_ana.csv'

def csv_dict_creator(input_csv, dict_out):
    input_csv = open(input_old_csv, 'r').read().splitlines()
    details_dict = {}
    for line in input_csv[1:]:
        col = line.split(',')
        cg_group = col[0].strip('"')
        details_dict[cg_group] = {}
        for cell in col[14:]:
            if '\t' not in cell and '""' not in cell:
                if 'ERR' in cell and "scaffold" not in cell:
                    group_n_gene = cell[:-4].strip('"')
                    strand_orein = cell[-3]
                    details_dict[cg_group][group_n_gene] = strand_orein
                elif "ERR" in cell and "scaffold" in cell:
                    group_n_gene = (cell.split('_')[0] + '_' + cell.split('_')[2][:-4]).strip('"')
                    strand_orein = cell[-3]
                    details_dict[cg_group][group_n_gene] = strand_orein
                elif "cds" in cell:
                    group_n_gene = "GCA000006745_" + cell.split('_')[0][4:].strip('"')
                    strand_orein = cell[-3]
                    details_dict[cg_group][group_n_gene] = strand_orein
                elif ".fna" in cell:
                    group_n_gene = "GCA" + cell.split('_')[1][:-2] + "_" + cell.split('_')[-1][:-4]
                    strand_orein = cell[-3]
                    details_dict[cg_group][group_n_gene] = strand_orein
                elif "GC" in cell:
                    group_n_gene = "GCA" + cell.split('_')[1] + "_" + cell.split('_')[-1][:-4]
                    strand_orein = cell[-3]
                    details_dict[cg_group][group_n_gene] = strand_orein

    # ##writing out dict for trouble shooting with ker errors
    # with open(dict_out, 'w') as file:
    #     json.dump(details_dict, file)
    return details_dict

def add_strand_orein_to_csv(input_csv, details_dict, outfile_path):
    input_csv = open(input_csv).read().splitlines()
    out_list = []

    for line in input_csv[1:]:
        new_line_list = []
        col = line.split(',')
        cg_group = col[0].strip('"')
        for cell in col:
            if '\t' not in cell and '_' in cell and "group" not in cell:
                group_n_gene = cell.split("_")[0] + "_" + cell.split("_")[2].strip("'")
                new_line_list.append(cell + '(' + details_dict[cg_group][group_n_gene] + ')')
            else:
                new_line_list.append(cell)
        out_list.append(new_line_list)

    with open(outfile_path, "w") as f:
        writer = csv.writer(f)
        writer.writerows(out_list)

# add_strand_orein_to_csv(input_cg_csv, csv_dict_creator(input_old_csv, dict_out_file), outfile_path)