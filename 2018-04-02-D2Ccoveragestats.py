import glob2
import numpy

outfile = open('/Users/liamcheneyy/Desktop/assembly_coverages11.csv', 'w')

for filename in glob2.iglob('/Users/liamcheneyy/Desktop/stat_files/*'):
    acc = filename.replace('/Users/liamcheneyy/Desktop/stat_files/','').replace('_assemblystats.txt','')
    infile = open(filename, 'r')
    outfile.write('\n')
    outfile.write(acc + '\n')
    avlen = []
    avcov = []
    # print(acc)
    for line in infile:
        if '>NODE' in line:
            split = line.split("_")
            len = int(split[3])
            avlen.append(len)
            cov = float(split[5])
            avcov.append(cov)
            avl = numpy.mean(avlen)
            avc = numpy.mean(avcov)
    print(str(avl) + '\t' + str(avc))
    # outfile.write('Average Length : ' + str(avl) + '\n')
    # outfile.write('Average Coverage : ' + str(avc) + '\n')

# outfile.close()
