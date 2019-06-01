import json
import os


PATHPARAM = 'CodeParameters/'
NKFILE = 'list_n_k.txt'
FILELIST = os.listdir(PATHPARAM)

print('Creating file: {}'.format(NKFILE))
with open(NKFILE, 'w') as nkfile:
    nkfile.write('n k\n')
    for FILE in FILELIST:
        print('Fetching and processing file: {}'.format(FILE))
        with open(PATHPARAM + FILE, 'r') as paramfile:
            jsondict = json.load(paramfile)
            nkfile.write('{} {}\n'.format(jsondict['n'], jsondict['k']))
