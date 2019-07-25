import json
import os


PATHPARAM = 'CodeParameters/'
FILELIST = [f for f in os.listdir(PATHPARAM)]

for FILE in FILELIST:
    with open(PATHPARAM + FILE, 'r') as paramfile:
        jsondict = json.load(paramfile)
        n, k = (jsondict['n'], jsondict['k'])
        if 'triorthogonal' in jsondict:
            trio = jsondict['triorthogonal']
            if trio and k > 1:
                print('Triorthogonal! file:{} -> {} {} {}\n'.format(FILE, n, k, jsondict["dz upper bound"]))
