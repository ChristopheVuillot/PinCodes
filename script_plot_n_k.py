import json
import os


PATHPARAM = 'CodeParameters/'
NKFILE = 'list_n_k_syst34_dim2.txt'
FILELIST = [f for f in os.listdir(PATHPARAM) if 'systematic34_dim2' in f]

print('Creating file: {}'.format(NKFILE))
with open(NKFILE, 'w') as nkfile:
    nkfile.write('n k d\n')
    for FILE in FILELIST:
        print('Fetching and processing file: {}'.format(FILE))
        with open(PATHPARAM + FILE, 'r') as paramfile:
            jsondict = json.load(paramfile)
            n, k = (jsondict['n'], jsondict['k'])
            if k > 0:
                nkfile.write('{} {} {}\n'.format(n, k, jsondict["dz upper bound"]))
            else:
                nkfile.write('{} {} {}\n'.format(n, k, 0))
