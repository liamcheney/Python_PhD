import sys

alleles_path = sys.argv[1]
project_path = sys.argv[2]
alleles_loc = '/Users/liamcheneyy/Alleles/Vibrio'
alleles_file = []
alleles_file = open(alleles_path).read().splitlines()

for num in range(0, len(alleles_file),1):
    if "alleleLocation = " in alleles_file[num]:
        alleles_file[num] = alleles_file[num].split('=')[0] + '= ' + '"' + alleles_loc + '"'

    if "projectpath = " in alleles_file[num]:
        alleles_file[num] = alleles_file[num].split('=')[0] + '= ' + '"' + project_path + '"'

    if 'addInfo(' in alleles_file[num]:
        alleles_file[num] = alleles_file[num].split('(')[0] + '(' + '"' + project_path + '"' + ',"Mgt",args.database,isolate_info)'

    if 'addTheHstMatrix(' in alleles_file[num]:
        alleles_file[num] = alleles_file[num].split('(')[0] + '(' + '"' +  project_path + '"' + ',"Mgt",args.database,mgtlist)'
    #
    if 'DbConString = ' in alleles_file[num]:
        alleles_file[num] = alleles_file[num].split('=')[0] + '= ' + '"' + "dbname='vcseventh' host='0.0.0.0' port='5432' user='postgres' password='Potle2368'" + '"' + '  ## connection info for Db - assign new user that can only do what is needed for script'

with open(alleles_path,'w') as out:
    for line in alleles_file:
        out.write(line + '\n')
