from time import sleep as sl
input = open('/Users/liamcheneyy/Desktop/ms_tree.nwk', 'r').read().splitlines()

acc_list = []
for i in input:
    everything = i.split(',')
    for j in everything:
        err = j.split(':')[0]
        if '(' in err:
            err = err.split('(')[-1]
        acc_list.append(err)

acc_list = set(acc_list)
#
for p in acc_list:
    print(p)