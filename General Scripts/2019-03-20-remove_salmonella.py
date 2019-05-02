import sys

settings_in_path = sys.argv[1]
mgt_path = sys.argv[2]

settings_in = []
settings_in = open(settings_in_path).read().splitlines()

for num in range(0, len(settings_in),1):
    if "    'Salmonella'," in settings_in[num]:
        settings_in[num] = ''
    if 'MEDIA_ROOT = "Uploads/"' in settings_in[num]:
        settings_in[num] = 'MEDIA_ROOT = ' + '"' + mgt_path + '"'
    if "'salmonella': {" in settings_in[num]:
        settings_in[num] = "    'vibrio': {"
    if "'NAME': 'salmonella50'," in settings_in[num]:
        settings_in[num] = "        'NAME': 'vcseventh',"
    if "APPS_DATABASE_MAPPING" in settings_in[num]:
        settings_in[num] = "APPS_DATABASE_MAPPING = {'Vibrio' : 'vibrio'}"
    if "                         'Vibrio': 'vibrio'}" in settings_in[num]:
        settings_in[num] = ''
    if "SUBDIR_ALLELES" in settings_in[num]:
        settings_in[num] = "SUBDIR_ALLELES = 'Alleles/Vibrio' # on testhgtdb & local"

del settings_in[35]

with open(settings_in_path,'w') as out:
    for line in settings_in:
        out.write(line + '\n')
