"""Script to analyse all codes
"""
from itertools import combinations
import os
import re
import json
import numpy as np
from QuantumCodeAnalysis.QuantumCodeAnalysis import low_weight_logical, logicals, distance_lower_bound
from QuantumCodeConstruction.utils import readsparsematrix
from flinalg import invert_permutation

PATHSYSTEMATICHP = 'PCMatrices/systematichp/'
PATH535 = 'PCMatrices/color_code_535/'
PATH353 = 'PCMatrices/color_code_353/'

RULE = re.compile(r'(.*)[XZ]\.sms')
FILESET = {RULE.match(f).group(1) for f in os.listdir(PATHSYSTEMATICHP)}

NTRIAL = 10

for FILEPREFIX in FILESET:

    PROPDICT = {}
    print('Computing properties of code ' + FILEPREFIX + ':')
    MX = readsparsematrix(PATHSYSTEMATICHP + FILEPREFIX + 'X.sms').todense()
    MZ = readsparsematrix(PATHSYSTEMATICHP + FILEPREFIX + 'Z.sms').todense()

    NX, NQ = MX.shape
    NZ, _ = MZ.shape

    CSSCOND = (np.dot(MX, MZ.transpose()) % 2).sum() == 0

    PROPDICT['CSS'] = CSSCOND
    PROPDICT['n'] = NQ
    PROPDICT['num X-checks'] = NX
    PROPDICT['num Z-checks'] = NZ

    XW = {MX[j, :].sum() for j in range(MX.shape[0])}
    ZW = {MZ[j, :].sum() for j in range(MZ.shape[0])}

    PROPDICT['X-checks weights'] = XW
    PROPDICT['Z-checks weights'] = ZW

    (LX, _), (LZ, _) = logicals(MX, MZ)
    if LX:
        LX = np.vstack(LX)
        LZ = np.vstack(LZ)
        GX = np.block([[LX], [MX]])
        GZ = np.block([[LZ], [MZ]])
        K, _ = LX.shape
        PROPDICT['k'] = K
    else:
        print('Zero logical operator !')
        PROPDICT['k'] = 0
        continue

    LOWX, PERMX = low_weight_logical(GX, LZ, NTRIAL)
    LOWZ, PERMZ = low_weight_logical(GZ, LX, NTRIAL)

    DXUP = (LOWX % 2).sum()
    DZUP = (LOWZ % 2).sum()

    PROPDICT['dx upper bound'] = DXUP
    PROPDICT['small X logical'] = LOWX[invert_permutation(PERMX)]
    PROPDICT['dz upper bound'] = DZUP
    PROPDICT['small Z logical'] = LOWZ[invert_permutation(PERMZ)]

    DXCOND = (np.dot(MZ[:, PERMX], LOWX) % 2).sum() == 0
    DXCOND = DXCOND and ((np.dot(LZ[:, PERMX], LOWX) % 2).sum() != 0)
    DZCOND = (np.dot(MX[:, PERMZ], LOWZ) % 2).sum() == 0
    DZCOND = DZCOND and ((np.dot(LX[:, PERMZ], LOWZ) % 2).sum() != 0)

    PROPDICT['valid dx upper bound'] = DXCOND
    PROPDICT['valid dz upper bound'] = DZCOND

    DZMORETHAN2, W2Z = distance_lower_bound(MX, LX, 2)
    DXMORETHAN2, W2X = distance_lower_bound(MZ, LZ, 2)

    PROPDICT['dz > 2'] = DZMORETHAN2
    if not DZMORETHAN2:
        PROPDICT['small Z logical'] = W2Z

    PROPDICT['dx > 2'] = DXMORETHAN2
    if not DXMORETHAN2:
        PROPDICT['small X logical'] = W2X

    TRICOND = True
    for ilog1, ilog2 in combinations(range(K), 2):
        for istab in range(NX):
            TEST = np.multiply(LX[ilog1, :],
                               np.multiply(LX[ilog2, :],
                                           MX[istab, :])).sum() % 2 == 0
            TRICOND = TRICOND and TEST
            if not TRICOND:
                break
        if not TRICOND:
            break
    PROPDICT['triorthogonal'] = TRICOND

    with open('CodeParameters/' + FILEPREFIX + '.json', 'w') as dictfile:
        dictfile.write(json.dump(PROPDICT))
