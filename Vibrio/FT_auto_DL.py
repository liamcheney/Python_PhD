import os
from time import sleep as sl
import shutil
import subprocess

def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-i", "--path", help="Path to main program directory. FUlL PATH REQUIRED.", default="/srv/scratch/lanlab/liam/2020-03-23-staph_sra/kdm_parallel/New/")
    parser.add_argument("-m", "--max_list", help="Maximum size of each subset list. MAX=2000", type=int, default=3)
    parser.add_argument("-x", "--max_kdm", help="Maximum number of KDM sessions.", type=int, default=5)
    parser.add_argument("-b", "--bash", help="Source liams bash profile.", default="/home/z5087966/.bash_profile")
    parser.add_argument("-f", "--force", help="Delete previous files.", default=True)
    parser.add_argument("-t", "--time", help="Time between checking for new KDM sessions (minutes).", default=.1)

    args = parser.parse_args()

    if args.max_list > 2000:
        raise argparse.ArgumentTypeError("Maximum subset size = 2000")
    return args

def downloader(accesion_lists_dict, complete_count, completed_list, print_help, args):
    active_windows = session_checker()
    if print_help == True:
        print("To monitor activate KDM sessions enter: watch screen -ls")
        print("To view last update on new KDM session view: " + args.path + "/Data/KDM_update.txt")
        print('')

    if (active_windows == 0) and (len(completed_list) > 1):
        print("Finished Downloading")
        exit()

    else:
        #start loop to continue while downloading
        for key in accesion_lists_dict.keys():
            if active_windows <= args.max_kdm: #check have not exceeded KDM screens number
                if key not in completed_list:
                    #create a screen under accession name and execute download script
                    print("Creating KDM session. Downloading subset : " + key)
                    proc = """screen -dmS {}; screen -S {}""".format(key,key) + " -X stuff '. " + args.path + "/Data/" + key + "/" + """prefetch.ssh\r'"""
                    execute = subprocess.Popen(proc,shell=True).wait()
                    sl(0.5)

                    completed_list.append(key)
                    active_windows = session_checker()
                    complete_count = complete_count + 1

            elif active_windows >= args.max_kdm:
                print_help = False
                with open(args.path + "/Data/KDM_update.txt", 'w') as out:

                    datetime_ob = datetime.datetime.now()
                    out.write("Checked at " + str(datetime_ob) + '\n')
                    out.write("After " + str(args.time) + " minutes.")

                secs = args.time * 60
                sl(secs)
                downloader(accesion_lists_dict, complete_count, completed_list, print_help, args)

def session_checker():
    proc = "screen -ls"
    process = subprocess.Popen(proc, shell=True, stdout=subprocess.PIPE)
    result = process.stdout.read()
    result = result.decode('utf8')
    res = result.split('\t')

    for i in res:
        if 'No Sockets' in i:
            num = 0
            return num
        elif ('Socket' in i) and ('No Sockets' not in i):
            num = int(i.split()[1])
            return num

def sub_dir_creator(accesion_lists_dict, args):

    for key, value in accesion_lists_dict.items():
        os.mkdir(args.path + "/Data/" + str(key))
        sub_list_out(key, value, args)
        ssh_command_out(key, args)

def ssh_command_out(key, args):

    shutil.copyfile(args.path + "/Data/scripts/prefetch.ssh", args.path + '/Data/' + str(key) + '/prefetch.ssh')

    proc1 = "sed -i 's/" + "REPLACE_KEY" + "/" + key + "/g' " + args.path + '/Data/' + str(key) + '/prefetch.ssh'
    execute1 = subprocess.Popen(proc1, shell=True)

def sub_list_out(key, value, args):
        with open(args.path + '/Data/' + key + '/' + str(key) + '_accessions.txt', 'w') as out:
            for i in value:
                out.write(i + '\n')

def file_spliter(args):

    infile = open(args.path + "/accesions.txt",'r').read().split()

    out_lists = {}
    infile_list_size = len(infile)
    number_of_list = int(infile_list_size/args.max_list) + 1

    list_count_low = 0
    list_count_high = args.max_list

    for i in range(0, number_of_list, 1):

        list_name = ""

        #set correct name for last list (not always an even of the max_list_size)
        if number_of_list - i != 1:
            list_name = str(list_count_low) + "_" + str(list_count_high)
            out_lists[list_name] = infile[list_count_low:list_count_high]

        elif number_of_list - i == 1:
            place_holder = len(str(args.max_list))
            num_split = [x for x in str(infile_list_size)]
            extra = num_split[-place_holder + 1:]
            extra_join = int(''.join(extra))
            list_count_high = list_count_high + extra_join
            list_name = str(list_count_low) + "_" + str(list_count_high)
            out_lists[list_name] = infile[list_count_low:list_count_high]

        #add for next out list
        list_count_low = list_count_low + args.max_list
        list_count_high = list_count_high + args.max_list


    keep_dict = {}
    for key, value in out_lists.items():
        if len(value) > 0:
            keep_dict[key] = value

    sorted_dict = dict(sorted(keep_dict.items()))

    return sorted_dict

def main():
    args = parseargs()
    print("Beginning Downloads. Sit back and relax =D")

    #check for exisiting files
    if args.force:
        import shutil
        shutil.rmtree(args.path + "/Data/")
    os.mkdir(args.path + "/Data/")

    #creates separate lists from input file
    accesion_lists_dict = file_spliter(args)

    #create sub-directories - adds ssh command
    sub_dir_creator(accesion_lists_dict, args)

    #downloading reads
    print("Maximum number of KDM sessions running : " + str(args.max_kdm))
    complete_count = 0
    completed_list = []
    print_help = True

    downloader(accesion_lists_dict, complete_count, completed_list, print_help, args)

if __name__ == '__main__':
    main()
