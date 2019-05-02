infile = open('/Users/liamcheneyy/Desktop/mlst_output.txt','r').read().splitlines()

st_dict = {}
for line in infile:
    if "vcholerae" in line:
        col = line.split("\t")
        ST = col[2]
        if ST not in st_dict.keys():
            st_dict[ST] = 1
        if ST in st_dict.keys():
            st_dict[ST] = st_dict[ST] + 1
        # print(col[0] + '\t' + col[2])

for key, value in st_dict.items():
    print(key + '\t' +  str(value))
