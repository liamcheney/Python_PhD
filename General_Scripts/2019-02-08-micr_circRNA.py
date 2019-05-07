infile = open('/Users/liamcheneyy/Desktop/circRNA_MIC(GW23).csv', 'r').read().splitlines()

mic_list = []
circRNA_list = []
for line in infile[1:145]:
    if ':-' not in line:
        col = line.split(',')
        mic_start = int(col[1].split(':')[-1].split('-')[0])
        mic_end = int(col[1].split(':')[-1].split('-')[-1])
        mic_name = col[0]

        list_n = [mic_start,mic_end, mic_name]
        mic_list.append(list_n)

for line in infile[1:19330]:
    col = line.split(',')
    circRNA_start = int(col[4])
    circRNA_end = int(col[5])
    circRNA_list.append([circRNA_start, circRNA_end])

yes_count = 0
no_count = 0
out_list = []
for i in mic_list:
    for j in circRNA_list:
        if i[0] > j[0] and i[1] < j[1]:
            yes_count = yes_count + 1
            out_list.append([i[2],i[0],i[1],j[0],j[1]])
        else:
            no_count = no_count + 1

outfile = open('/Users/liamcheneyy/Desktop/results_circRNA_MIC(GW23).csv', 'w')
outfile.write('MicroExon' + ',' + 'Micro Start' + ',' + 'Micro End' + ',' + 'CircRNA Start' + ',' + 'CircRNA End' + '\n')
for j in out_list:
    for k in j:
        outfile.write(str(k) + ',')
    outfile.write('\n')
outfile.close()