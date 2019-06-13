"""Script to analyse all codes
"""
from itertools import combinations
import os
import re
import json
import numpy as np
from QuantumCodeAnalysis.QuantumCodeAnalysis import low_weight_logical, logicals  # , distance_lower_bound
from QuantumCodeConstruction.utils import readsparsematrix
from flinalg import invert_permutation


PATHSYSTEMATICHP = 'PCMatrices/systematichp/'
PATHSYSTEMATICHP43 = 'PCMatrices/systematichp43/toanalyse/'
PATHNARROWCC = 'PCMatrices/narrowCC/'
PATH535 = 'PCMatrices/color_code_535/'
PATH353 = 'PCMatrices/color_code_353/'

PATH = PATHSYSTEMATICHP43

RULE = re.compile(r'(.*)[XZ]\.sms')
FILESET = {RULE.match(f).group(1) for f in os.listdir(PATH) if '(33)' in f}

NTRIAL = 30

for FILEPREFIX in FILESET:

    PROPDICT = {}
    print('Computing properties of code ' + FILEPREFIX + ':')
    print('Loading check matrices')
    MX = readsparsematrix(PATH + FILEPREFIX + 'X.sms').todense()
    MZ = readsparsematrix(PATH + FILEPREFIX + 'Z.sms').todense()

    NX, NQ = MX.shape
    NZ, _ = MZ.shape

    print('Checking CSS condition')
    CSSCOND = (np.dot(MX, MZ.transpose()) % 2).sum() == 0

    PROPDICT['CSS'] = int(CSSCOND)
    PROPDICT['n'] = int(NQ)
    PROPDICT['num X-checks'] = int(NX)
    PROPDICT['num Z-checks'] = int(NZ)

    print('computing stabilizer weigths')
    XW = {MX[j, :].sum() for j in range(MX.shape[0])}
    ZW = {MZ[j, :].sum() for j in range(MZ.shape[0])}

    PROPDICT['X-checks weights'] = [int(w) for w in XW]
    PROPDICT['Z-checks weights'] = [int(w) for w in ZW]

    print('Computing logical operators')
    (LX, _), (LZ, _) = logicals(MX, MZ)
    if LX:
        LX = np.vstack(LX)
        LZ = np.vstack(LZ)
        GX = np.block([[LX], [MX]])
        GZ = np.block([[LZ], [MZ]])
        K, _ = LX.shape
        PROPDICT['k'] = int(K)
    else:
        print('Zero logical operator !')
        PROPDICT['k'] = int(0)
        with open('CodeParameters/' + FILEPREFIX + '.txt', 'w') as dictfile:
            print('writing properties to file')
            dictfile.write(json.dumps(PROPDICT))
        continue

    print('Finding distance upper bound')
    LOWX, PERMX = low_weight_logical(GX, LZ, NTRIAL)
    LOWZ, PERMZ = low_weight_logical(GZ, LX, NTRIAL)

    DXUP = (LOWX % 2).sum()
    DZUP = (LOWZ % 2).sum()

    PROPDICT['dx upper bound'] = int(DXUP)
    PROPDICT['small X logical'] = [int(b) for b in LOWX[invert_permutation(PERMX)]]
    PROPDICT['dz upper bound'] = int(DZUP)
    PROPDICT['small Z logical'] = [int(b) for b in LOWZ[invert_permutation(PERMZ)]]

    DXCOND = (np.dot(MZ[:, PERMX], LOWX) % 2).sum() == 0
    DXCOND = DXCOND and ((np.dot(LZ[:, PERMX], LOWX) % 2).sum() != 0)
    DZCOND = (np.dot(MX[:, PERMZ], LOWZ) % 2).sum() == 0
    DZCOND = DZCOND and ((np.dot(LX[:, PERMZ], LOWZ) % 2).sum() != 0)

    PROPDICT['valid dx upper bound'] = int(DXCOND)
    PROPDICT['valid dz upper bound'] = int(DZCOND)

    # print('checking distance larger than 2')
    # DZMORETHAN2, W2Z = distance_lower_bound(MX, LX, 2)
    # DXMORETHAN2, W2X = distance_lower_bound(MZ, LZ, 2)

    # PROPDICT['dz > 2'] = int(DZMORETHAN2)
    # if not DZMORETHAN2:
    #     PROPDICT['small Z logical'] = list(W2Z)

    # PROPDICT['dx > 2'] = int(DXMORETHAN2)
    # if not DXMORETHAN2:
    #     PROPDICT['small X logical'] = list(W2X)

    # print('Checking triorthogonality')
    # TRICOND = True
    # for ilog1, ilog2 in combinations(range(K), 2):
    #     for istab in range(NX):
    #         TEST = np.multiply(LX[ilog1, :],
    #                            np.multiply(LX[ilog2, :],
    #                                        MX[istab, :])).sum() % 2 == 0
    #         TRICOND = TRICOND and TEST
    #         if not TRICOND:
    #             break
    #     if not TRICOND:
    #         break
    # PROPDICT['triorthogonal'] = int(TRICOND)

    with open('CodeParameters/' + FILEPREFIX + '.txt', 'w') as dictfile:
        print('writing properties to file')
        dictfile.write(json.dumps(PROPDICT))
