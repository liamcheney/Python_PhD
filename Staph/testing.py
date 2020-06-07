import argparse
from time import sleep as sl
import multiprocessing as mp


def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("-d", "--database_name", required=True,
    #                     help="sql database to search (eg. vcseventh)")

    args = parser.parse_args()

    return args

def process(line):

    pass

def main():
    args = parseargs()

    # init objects
    pool = mp.Pool(mp.cpu_count())
    jobs = []

    # create jobs
    with open("/Users/liamcheneyy/Desktop/out.txt") as f:
        for line in f:
            jobs.append(pool.apply_async(process(line), (line)))

    # wait for all jobs to finish
    for job in jobs:
        job.get()

    # clean up
    pool.close()

if __name__ == '__main__':
    main()
