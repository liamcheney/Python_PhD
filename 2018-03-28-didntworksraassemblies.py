import glob2

lowcov = []
lowmem = []
for filename in glob2.iglob('/Users/liam/Desktop/didnt_work/*'):
    fl = filename.replace('/Users/liam/Desktop/didnt_work/', '').replace('_spades.log', ' ')
    with open(filename, 'r') as infile:
        for line in infile:
            if 'he reads contain too many k-mers' in line:
                print(fl)
            if 'make sure that the coverage is indeed uniform' in line:
                print(fl)
