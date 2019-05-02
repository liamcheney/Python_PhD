import glob
from time import sleep as sl
import sys

# glob_sum_in = sys.argv[1]
glob_sum_in = "/Users/liamcheneyy/Desktop/8_Snap_extra/*"

for filename in glob.iglob(glob_sum_in):
    name = filename.split('/')[-1].strip('summary.')
    with open(filename, 'r') as file:
        for line in file:
            if "Averages of all pairwise comparisons" in line:
                ds = float(line.split(':')[1].split(',')[0].split('=')[-1].strip('  '))
                dn = float(line.split(':')[1].split(',')[1].split('=')[-1].strip('  '))
                dn_ds = round(1 / float(line.split(':')[1].split(',')[-2].split('=')[-1].strip('  ')),8)

                print(name, dn, ds, dn_ds)
