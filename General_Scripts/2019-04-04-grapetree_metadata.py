from time import sleep as sl

#infiles
mgt_st_path = "/Users/liamcheneyy/Desktop/MGT_sts.tsv"
metadata_path = "/Users/liamcheneyy/MGT_files/metadata/processed_metadata9.tsv"
meta_want_list = ["Year","Country","Continent","Wave","Serogroup","Biotype","Serotype","Strain","UN Subregion","City","State"]
outpath = "/Users/liamcheneyy/Desktop/Grapetree_metadata.txt"

def read_in_make_dicts(mgt_st_path, metadata_path):
    print("Reading in MGT STs file.")

    mgt_st_dict = {}
    mgt_st_file = open(mgt_st_path,'r').read().splitlines()
    for line in range(1, len(mgt_st_file), 1):
        col = mgt_st_file[line].split('\t')
        run = col[0].strip()
        mgt_st_dict[run] = {}
        for cell in range(0, len(col), 1):
            for category_line in range(0, 1, 1):
                header = mgt_st_file[category_line].split('\t')
                mgt_st_dict[run][header[cell]] = col[cell].upper()

    print("Reading in Metadata File.")

    meta_data_dict = {}
    meta_file = open(metadata_path, 'r').read().splitlines()
    for line in range(1, len(meta_file), 1):
        col = meta_file[line].split('\t')
        run = col[0].strip()
        meta_data_dict[run] = {}
        for cell in range(0, len(col), 1):
            for category_line in range(0, 1, 1):
                header = meta_file[category_line].split('\t')
                meta_data_dict[run][header[cell]] = col[cell].upper()

    return mgt_st_dict, meta_data_dict

def clean_up(info_dict, new_file_dict):
    #will go over dictionaries and remove any cells not providing metadata such as "not assigned".
    for key, value in info_dict.items():
        for item in info_dict[key]:
            if "UNKNOWN" in info_dict[key][item]:
                info_dict[key][item] = ''
            if "NOT" in info_dict[key][item]:
                info_dict[key][item] = ''
            if "NONE" in info_dict[key][item]:
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
            if "NONE" in new_file_dict[key][item]:
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

def choose_metadata(mgt_st_dict, meta_data_dict, meta_want_list):
    print("Finding and adding desired metadata.")
    final_metadata_dict={}
    no_meta_count = 0
    for key in mgt_st_dict.keys():
        for element in meta_want_list:
            try:
                mgt_st_dict[key][element] = meta_data_dict[key][element]
            except:
                no_meta_count = no_meta_count + 1
                continue

    print("No metadata could be found for " + str(no_meta_count) + " strains.")
    final_metadata_dict = mgt_st_dict

    return final_metadata_dict

def write_out(final_metadata_dict,out_path, meta_want_list):
    print("Writing out data.")
    with open(out_path,'w') as out:
        #write out file header
        key_h = list(final_metadata_dict.keys())[0]

        for item in final_metadata_dict[key_h]:
            out.write(item + '\t')
        for i in meta_want_list:
            out.write(i + '\t')
        out.write('\n')

        #write out content
        for key, value in final_metadata_dict.items():
            for cell in final_metadata_dict[key]:
                out.write(final_metadata_dict[key][cell] + '\t')
            out.write('\n')

def main(mgt_st_path, meta_want_list, outpath):

    #create dicts from inputs
    mgt_st_dict, meta_data_dict = read_in_make_dicts(mgt_st_path, metadata_path)

    #clean up
    mgt_st_dict, meta_data_dict = clean_up(mgt_st_dict, meta_data_dict)

    #choose the metadata (in addition to all MGT STs) to be added
    final_metadata_dict = choose_metadata(mgt_st_dict, meta_data_dict, meta_want_list)

    #write out
    write_out(final_metadata_dict, outpath, meta_want_list)

main(mgt_st_path, meta_want_list, outpath)