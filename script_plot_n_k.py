import json
import os


PATHPARAM = 'CodeParameters/'
NKFILE = 'list_n_k_narrowCC_33_dim6.txt'
FILELIST = [f for f in os.listdir(PATHPARAM) if 'narrowCC' in f and 'dim6' in f and '(33)' in f]

print('Creating file: {}'.format(NKFILE))
with open(NKFILE, 'w') as nkfile:
    nkfile.write('n k d wx wz\n')
    for FILE in FILELIST:
        print('Fetching and processing file: {}'.format(FILE))
        with open(PATHPARAM + FILE, 'r') as paramfile:
            jsondict = json.load(paramfile)
            n, k = (jsondict['n'], jsondict['k'])
            wx, wz = max(jsondict["X-checks weights"]), max(jsondict["Z-checks weights"])
            if k > 0:
                nkfile.write('{} {} {} {} {}\n'.format(n, k, jsondict["dz upper bound"], wx, wz))
            else:
                nkfile.write('{} {} {} {} {}\n'.format(n, k, 0, wx, wz))
