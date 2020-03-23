from time import sleep as sl
import sys
import glob

#input files
infile = "/Users/liamcheneyy/MGT_files/metadata/processed_metadata9.tsv"
new_file = "/Users/liamcheneyy/Desktop/metadata.txt"
out_path = "/Users/liamcheneyy/Desktop/metadata/processed_metadata9.tsv"
countries = "/Users/liamcheneyy/MGT_files/metadata/country-to-cont.tsv"

def read_in(infile, new_file):

    infile = open(infile,'r').read().splitlines()
    new_file = open(new_file,'r').read().splitlines()

    return infile, new_file

def create_country_to_continent_dict(countries):

    inlist = open(countries).read().splitlines()
    countries_dict = {}

    for line in inlist:
        col = line.split('\t')
        countries_dict[col[1].strip('\xa0')] = col[0].strip('\xa0')
    return countries_dict

def create_dicts_from_files(infile, new_file):

    print("Reading in all metadata files.")
    ###create template file to add all new data to
    info_dict = {}
    for line in range(1, len(infile), 1):
        col = infile[line].split('\t')
        run = col[0].strip()
        info_dict[run] = {}
        for cell in range(0, len(col), 1):
            for category_line in range(0, 1, 1):
                header = infile[category_line].split('\t')
                info_dict[run][header[cell]] = col[cell].upper()

    ###create dicts for the remaining files
    new_file_dict = {}
    for line in range(1, len(new_file), 1):
        col = new_file[line].split('\t')
        run = col[0].strip()
        new_file_dict[run] = {}
        for cell in range(0, len(col), 1):
            for category_line in range(0, 1, 1):
                header = new_file[category_line].split('\t')
                new_file_dict[run][header[cell]] = col[cell].upper()

    return info_dict, new_file_dict

def clean_up(info_dict, new_file_dict):
    #will go over dictionaries and remove any cells not providing metadata such as "not assigned".
    for key, value in info_dict.items():
        for item in info_dict[key]:
            if "UNKNOWN" in info_dict[key][item]:
                info_dict[key][item] = ''
            if "NOT" in info_dict[key][item]:
                info_dict[key][item] = ''
            if "N/A" in info_dict[key][item]:
                info_dict[key][item] = ''
            if "NAN" in info_dict[key][item]:
                info_dict[key][item] = ''
            if "MISSING" in info_dict[key][item]:
                info_dict[key][item] = ''
            if "/" in info_dict[key][item]:
                info_dict[key][item] = ''
            if info_dict[key][item] == "WT":
                info_dict[key][item] = 'WILD TYPE'
            if "CLINICAL" in info_dict[key][item]:
                info_dict[key][item] = 'HUMAN'
            if "NO D" in info_dict[key][item]:
                info_dict[key][item] = ''
            if info_dict[key][item] == "ENVIRONMENTAL":
                info_dict[key][item] = 'WATER'
            if "VIET" in info_dict[key][item]:
                info_dict[key][item] = 'VIETNAM'
            if "REPUBLIC OF SOUTH AFRICA" in info_dict[key][item]:
                info_dict[key][item] = 'SOUTH AFRICA'

    for key, value in new_file_dict.items():
        for item in new_file_dict[key]:
            if "UNKNOWN" in new_file_dict[key][item]:
                new_file_dict[key][item] = ''
            if "NOT" in new_file_dict[key][item]:
                new_file_dict[key][item] = ''
            if "N/A" in new_file_dict[key][item]:
                new_file_dict[key][item] = ''
            if "NAN" in new_file_dict[key][item]:
                new_file_dict[key][item] = ''
            if "MISSING" in new_file_dict[key][item]:
                new_file_dict[key][item] = ''
            if "/" in new_file_dict[key][item]:
                new_file_dict[key][item] = ''
            if "CLINICAL" in new_file_dict[key][item]:
                new_file_dict[key][item] = 'HUMAN'
            if "NO D" in new_file_dict[key][item]:
                new_file_dict[key][item] = ''
            if new_file_dict[key][item] == "WT":
                new_file_dict[key][item] = 'WILD TYPE'
            if new_file_dict[key][item] == "ENVIRONMENTAL":
                new_file_dict[key][item] = 'WATER'
            if "VIET" in new_file_dict[key][item]:
                new_file_dict[key][item] = 'VIETNAM'
            if "REPUBLIC OF SOUTH AFRICA" in new_file_dict[key][item]:
                new_file_dict[key][item] = 'SOUTH AFRICA'




    return info_dict, new_file_dict

def add_new_strains(info_dict, new_file_dict, countries, ):

    print("Adding new data.")
    ##if the strain does not already exist
    info_dict, new_file_dict = if_strain_does_not_exist(info_dict, new_file_dict)

    ##if a strain already exists in the dictionary
    info_dict = if_strain_exists(info_dict, new_file_dict, countries)
    return info_dict

def check_new_categories(info_dict, infile, new_file_dict, new_file):

    #make list of existing categories from template and newfile
    existing_cat_list = [x for x in infile[0].split('\t')]

    #make list of categories from new file
    new_file_cats = [x for x in new_file[0].split('\t')]

    #check if new cat is not in existing cats
    for cat in new_file_cats:
        if cat not in existing_cat_list:
            found_new_cat = 2
            # found_new_cat = int(input(cat + '\t' + '\t' + '\t' + "New metadata cateogry found." + '\t' + '\t' + '\t' + "Return 1: to skip adding this category of metadata, Return 2: to add new category, Return 3: to change the new files category name."))

            ##skip new found category
            if found_new_cat == 1:
                for key, value in new_file_dict.items():
                    del new_file_dict[key][cat]

            ##add new category to dictionary
            if found_new_cat == 2:
                for key, value in info_dict.items():
                    info_dict[key][cat] = ''
                    existing_cat_list.append(cat)

            #edit an existing metadata category
            if found_new_cat == 3:

                #check if new label is in dict
                key_to_change = str(input("Existing = " + str(cat) + '\t' + ". Change new category label to what ... ?"))

                #check if new label is apart of template
                if key_to_change not in existing_cat_list:
                    print('\t' + '\t' + '\t' + "ERROR:" + '\t' + "Please rename label to existing category in template file.")
                    print('\t' + '\t' + '\t' + "ERROR: " + str(cat) + '\t' + " not found in tempalte file.")
                    key_to_change = str(input('\t' + '\t' + '\t' + "Have another go : "))

                for key, value in new_file_dict.items():
                    new_file_dict[key][key_to_change] = new_file_dict[key].pop(cat)

    return info_dict, new_file_dict, existing_cat_list

def if_strain_exists(info_dict, new_file_dict, countries):
    country_dict = create_country_to_continent_dict(countries)

    ##create list of categories in new_file
    new_file_keys = []
    for i in new_file_dict.keys():
        for x in new_file_dict[i]:
            new_file_keys.append(x)

    print("New data found for accesions :")
    for key, value in info_dict.items():
        #checking the strain is in both template and new_file

        if key in new_file_dict.keys():
            for item in info_dict[key]:

                #if the template is missing data for a category
                if info_dict[key][item] == '':

                    ##find the data in the newfile
                    if item in new_file_keys:

                        ##if the new file is missing then pass
                        if new_file_dict[key][item] == '':
                            pass

                        ##if new file has data add to template
                        if new_file_dict[key][item] != '':
                            info_dict[key][item] = new_file_dict[key][item]

        #add continent automatically
        try:
            if info_dict[key]['Country'] != '':
                info_dict[key]['Continent'] = country_dict[info_dict[key]['Country']]
        except:
            continue

    #         ##if the template already has data
    #         if info_dict[key][item] != '':
    #
    #             # check the new file has the datatype
    #             if item in new_file_keys:
    #
    #                 # if the new data is different to the existing data
    #                 # if new_file_dict[key][item].strip() != info_dict[key][item].strip():
    #                 #     conflicting_data_dict[key][item] = [info_dict[key][item].strip(), new_file_dict[key][item].strip()]
    #
    #                 # if template and new are the same then skip
    #                 if new_file_dict[key][item].strip() == info_dict[key][item].strip():
    #                     pass



    return info_dict

def if_strain_does_not_exist(info_dict, new_file_dict):
    #some strains have genbank IDs instead

    #identify strain not found in templae_dict
    for key in new_file_dict.keys():
        if key not in info_dict.keys():

            #some strains use a SRA_experiment/accession instead, convert to Run ID
            for temp_key in info_dict.keys():
                if info_dict[temp_key]['SRA_Sample'] == key:
                    new_file_dict[temp_key] = new_file_dict.pop(key)

                try:
                    if info_dict[temp_key]['SRA_accession'] == key:
                        new_file_dict[temp_key] = new_file_dict.pop(key)
                except:
                    continue

    return info_dict, new_file_dict

def write_out(out_path, info_dict):
    print("Writing out data.")
    with open(out_path,'w') as out:
        #write out file header
        key_h = list(info_dict.keys())[0]
        for item in info_dict[key_h]:
            out.write(item + '\t')
        out.write('\n')

        #write out content
        for key, value in info_dict.items():
            # print(key)
            for j in info_dict[key]:
                if isinstance(info_dict[key][j], (list,)):
                    comb = ','.join(info_dict[key][j])
                    out.write(comb + '\t')
                else:
                    out.write(info_dict[key][j] + '\t')
            out.write('\n')

    ##writout the conflicting dict
    # with open(out_path, 'w') as out:

def main(infile, new_file, out_path, countries):

    ##read in files
    infile, new_file = read_in(infile, new_file)

    ##create template dict from SRA template data and new file
    info_dict, new_file_dict = create_dicts_from_files(infile, new_file)

    ##clean up both dicts to remove "not assigned, unknown etc"
    info_dict, new_file_dict = clean_up(info_dict, new_file_dict)

    ##search new file to find and add new categories
    info_dict, new_file_dict, existing_cat_list = check_new_categories(info_dict, infile, new_file_dict, new_file)

    ##will iterate over new information and add to info_dict{}.
    info_dict = add_new_strains(info_dict, new_file_dict, countries)

    #write out
    write_out(out_path, info_dict)

main(infile, new_file, out_path, countries)