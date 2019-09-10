import argparse
from time import sleep as sl

def read_in(infile_path):

    infile = open(infile_path, 'r').read().splitlines()
    return infile

def processing(infile):

    regions_list = []

    for line in infile[1:]:
        start_new_region = False

        col = line.split(',')

        #keep the first start
        first_start = 1
        #handle hits on reverse strand
        current_start, current_end = rev_orientation(col)

        #finding overlapping regions
        for line2 in infile[1:]:
            col2 = line2.split(',')
            new_start, new_end = rev_orientation(col2)
            start_new_region = False

            if start_new_region == True:
                break

            #check if next region is overlapping
            elif new_start < current_end:

                #if new frag has larger end, then extend
                if new_end > current_end:
                    current_end = new_end

            elif new_start > current_end:
                regions_list.append([current_start, current_end])
                # print(current_start, current_end)
                start_new_region = True


        # print(current_start , current_end)



    print(regions_list)



def rev_orientation(col):
    if int(col[1]) > int(col[2]):
        start = int(col[2])
        end = int(col[1])

    else:
        start = int(col[1])
        end = int(col[2])

    return start,end


def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args


def main(infile_path):
    args = parseargs()

    #read in data file
    infile = read_in(infile_path)

    #process infile
    output_list = processing(infile)

    #write oute


if __name__ == '__main__':

    infile_path = '/Users/liamcheneyy/Desktop/Book1.csv'
    main(infile_path)
