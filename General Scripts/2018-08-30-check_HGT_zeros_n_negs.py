from time import sleep as sl
import glob
from openpyxl import  Workbook


##go through the gene_profile for st
##save al the profiles for that st
##do this for each scheme

st_gene_profiles_folder = '/Users/liam/Desktop/Vibrio_D3A_hierarchical_ST/'
gene_anno_folder = '/Users/liam/Desktop/Vibrio_D3A_hierarchical_ST/'
scheme_selection = reversed(["MGT3"]) ##change based on scheme wanted

def check_st_data(st_gene_profiles_folder, gene_anno_folder, scheme_selection):
    results_dict = {}
    wb = Workbook()
    ws = wb.active

    for filename in glob.iglob(st_gene_profiles_folder + "*_gene_profiles.txt"):
        infile = open(filename, 'r').read().splitlines()
        scheme_num = filename.split('/')[-1].split('_')[0]
        results_dict[scheme_num] = {}

        #create dict with key = ST, value = profile
        st_prof  = {}
        for line in infile[1:]:
            col = line.split('\t')
            st = col[0]
            st_prof[st] = list(col[2:])
            st_zero_count = st_prof[st].count("0")

            st_neg_count = 0
            for x in st_prof[st]:
                if "-" in x:
                    st_neg_count = st_neg_count + 1
            results_dict[scheme_num][st] = {"ST":st, "Zero":st_zero_count, "Neg":st_neg_count}


    for scheme in scheme_selection:
        ws = wb.create_sheet(scheme)
        print(scheme)
        print("ST" + '\t' + "Zero Count" + '\t' + "Neg Count")
        for seq_typ in results_dict[scheme].keys():
            print(seq_typ + '\t' + str(results_dict[scheme][seq_typ]["Zero"]) + '\t' + str(results_dict[scheme][seq_typ]["Neg"]))
        print('\n')

    wb.save('/Users/liam/Desktop/test.xlsx')



check_st_data(st_gene_profiles_folder, gene_anno_folder,scheme_selection)
