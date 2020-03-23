import numpy
import tifffile
import oiffile
import pty

image = '/Users/liamcheneyy/Desktop/24h_1191_1A_60x.oib'
oib_path = '/Users/liamcheneyy/Desktop/oibs/'

try:
    oiffile.oib2oif(image, location=oib_path)
except:
    print("metadata file already created.")


infile = open('/Users/liamcheneyy/Desktop/oibs/24h_1191_1_60x.oif.files/s_C001Z001.pty').read().splitlines()
print(infile)

