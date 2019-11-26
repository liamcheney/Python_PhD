import pandas as pd
import sys
from time import sleep as sl

def group_genes(matrix_input,gene_function_dict,output):
    dic_pre_ab = matrix_input.T
    col = dic_pre_ab.iloc[0,:].tolist()
    dic_pre_ab.columns = col
    dic_pre_ab = dic_pre_ab.iloc[1:,:]
    dic_pre_ab = dic_pre_ab.to_dict()

    pre_ab_dic1 = {}


    for x, y in dic_pre_ab.items():
        pre_ab_dic1[x] = {}
        for a,b in y.items():
            if b>0:
                pre_ab_dic1[x][a] = b
    pre_ab_dic2={}
    for k, dic1 in pre_ab_dic1.items():
        pre_ab_dic2[k] = {}
        for gene, e in dic1.items():
            gene_fun = gene_function_dict[gene]
            pre_ab_dic2[k][gene_fun] = []
    for k, dic1 in pre_ab_dic1.items():
        for gene, e in dic1.items():
            gene_fun = gene_function_dict[gene]
            pre_ab_dic2[k][gene_fun].append(gene)
    ### pre_ab_dic1 refers to the genes presented, pre_ab_dic2 refers to the resistant functions harbored
    num_gene_acc = {}
    for c, d in pre_ab_dic1.items():
        num_gene_acc[c] = len(d)
    num_fun_acc = {}
    for k, v in pre_ab_dic2.items():
        num_fun_acc[k] = len(v)
    df_summary = pd.DataFrame.from_dict(num_gene_acc, orient='index')
    df_summary.columns = ['num_gene']
    df_summary['num_AR_type'] = df_summary.index.map(num_fun_acc)
    df_summary['num_type_detail'] = df_summary.index.map(pre_ab_dic2)
    df_summary.to_csv(output)
    print(df_summary.iloc[0:5,0:5])
def create_gene_dict(card_db_path,matrix_input):

    #create dict of card database
    card_db_dict = {}
    card_db = open(card_db_path, 'r').read().splitlines()

    #create dict for classes
    class_dict = {}

    for line in card_db[1:]:
        col = line.split(',')
        Gene_Name = col[0]
        DNA_Accession = col[1]
        AMR_Gene_Family = col[2]
        Drug_Class = col[3]
        group = col[4]
        Resistance_Mechanism = col[5]
        MDRC = col[6]

        card_db_dict[Gene_Name] = {'dna_acc':DNA_Accession,'AMR_fam':AMR_Gene_Family,'drug_class':Drug_Class,
                                   'resis_mech':Resistance_Mechanism, 'MDRC':MDRC, 'group':group}

        # class_dict[]
    #read in matrix
    matrix = pd.read_csv(matrix_input, index_col=0)

    return card_db_dict,matrix,class_dict
def calc_genes_per_strains(card_db_dict,matrix):

    gene_info_dict = {}
    AMR_list = list(matrix.columns.values)

    #everything per strain
    for index, row in matrix.iterrows():
        gene_info_dict[index] = {}

        #get list of AMR genes
        gene_list = list(row[row > 0].index)
        gene_info_dict[index]['genes'] = gene_list

        #get number of genes
        gene_count = row[row > 0].shape[0]
        gene_info_dict[index]['gene_count'] = gene_count

        #get number of drug classes etc
        drug_class_list = []
        AMR_fam_list = []
        resis_mech_list = []
        group_list = []

        for gene in gene_list:
            drug_class = card_db_dict[gene]['drug_class']
            resis_mech = card_db_dict[gene]['resis_mech']
            AMR_fam = card_db_dict[gene]['AMR_fam']
            group = card_db_dict[gene]['group']

            if drug_class not in drug_class_list:
                drug_class_list.append(drug_class)

            if AMR_fam not in AMR_fam_list:
                AMR_fam_list.append(AMR_fam)

            if resis_mech not in resis_mech_list:
                resis_mech_list.append(resis_mech)

            if group not in group_list:
                group_list.append(group)

        print(index, gene_count, len(drug_class_list), len(AMR_fam_list), len(resis_mech_list), len(group_list),  sep=',')

        gene_info_dict[index]['class_count'] = len(drug_class_list)
        gene_info_dict[index]['fam_count'] = len(AMR_fam_list)
        gene_info_dict[index]['resis_count'] = len(resis_mech_list)

        # print(gene_info_dict[index])
        # sl(1)







    # for line in matrix[1:2]:
    #     col = line.split(',')
    #     strain = col[0]
    #     gene_info_dict[strain] = {'gene_list':[]}
    #     gene_fam_count = []
    #
    #     gene_num = len([x for x in col[1:] if int(x) > 0])
    #
    #     for gene in range(1, len(col[1:]) + 1,1):
    #
    #         if int(line[gene]) > 0:
    #             gene_info_dict[strain]['gene_list'].append(line[gene])
    #
    # print(gene_info_dict)




    # number of other classes


def main():
    #paths
    matrix_input = '/Users/liamcheneyy/Desktop/amr/AMR_matrix.csv'
    output_path = '/Users/liamcheneyy/Desktop/amr/test_AMR_matrix.csv'
    card_db_path = '/Users/liamcheneyy/Desktop/amr/card_database.csv'

    #create gene functions dict
    card_db_dict,matrix,class_dict = create_gene_dict(card_db_path,matrix_input)

    #count AMR per strains
    save_dict = calc_genes_per_strains(card_db_dict,matrix)

    #output

    #TODO get list for all classes

if __name__ == '__main__':
    main()

