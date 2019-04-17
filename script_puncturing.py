"""script for puncturing codes
"""
# import gc
import numpy as np
import flinalg as fl
from QuantumCodeAnalysis.QuantumCodeAnalysis import low_weight_logical, distance_lower_bound
from QuantumCodeAnalysis.puncturing import puncture
from QuantumCodeConstruction.utils import readsparsematrix
# from permutations import PERM2222244, K2222244


MX = readsparsematrix('PCMatrices/systematichp/systematic22_dim6_transpose_1_X.sms').todense()
MZ = readsparsematrix('PCMatrices/systematichp/systematic22_dim6_transpose_1_Z.sms').todense()
# MX = readsparsematrix('PCMatrices/narrowCC/narrowCC2_dim6_X.sms').todense()
# MZ = readsparsematrix('PCMatrices/narrowCC/narrowCC2_dim6_Z.sms').todense()
# MX = readsparsematrix('PCMatrices/narrowCC/narrowCC_2222244_dim6_X.sms').todense()
# MZ = readsparsematrix('PCMatrices/narrowCC/narrowCC_2222244_dim6_Z.sms').todense()
# MX = readsparsematrix('PCMatrices/535_3420_XCOS.sms').todense()
# MZ = readsparsematrix('PCMatrices/535_3420_ZCOS.sms').todense()


RX, NQ = MX.shape
RZ, _ = MZ.shape

# MX = np.array(MX, dtype='uint8')
# SX = fl.row_reduce_transform(MX)

BESTGAMMA = 5
for k in range(100, 301):
    print('trying {} punctures:'.format(k), flush=True)
    for _ in range(1):
        PERM = np.random.permutation(NQ)
        # PERM = PERM2222244
        K = k

        PUNCTMATX, PERMSTD = puncture(MX, PERM, K)

        PUNCTSTABX = np.array(PUNCTMATX[K:, :], dtype='uint8')
        PUNCTLOGX = np.array(PUNCTMATX[:K, :], dtype='uint8')

        PUNCTSTABZ = fl.kernel(PUNCTMATX.transpose())
        KERZ = fl.kernel(PUNCTSTABX.transpose())
        # KERX = fl.kernel(PUNCTSTABZ.transpose())

        PUNCTLOGZ, _ = fl.quotient_basis(KERZ, PUNCTSTABZ)
        # REALLOGX, _ = fl.quotient_basis(KERX, PUNCTSTABX)

        TRIALS = 30
        LOWWEIGHTLOGZ, _ = low_weight_logical(KERZ, PUNCTLOGX, TRIALS)
        WEIGHT = LOWWEIGHTLOGZ.sum()
        KP = len(PUNCTLOGZ)
        NP = PUNCTMATX.shape[1]
        if WEIGHT > 1:
            GAMMA = np.log(NP / KP) / np.log(WEIGHT)
        else:
            GAMMA = 10
        if GAMMA < BESTGAMMA:
            TRIALS = 100
            LOWWEIGHTLOGZ, _ = low_weight_logical(KERZ, PUNCTLOGX, TRIALS)
            if WEIGHT > LOWWEIGHTLOGZ.sum():
                WEIGHT = LOWWEIGHTLOGZ.sum()
            if WEIGHT > 1:
                GAMMA = np.log(NP / KP) / np.log(WEIGHT)
            else:
                GAMMA = 10

        if GAMMA < BESTGAMMA:
            BESTGAMMA = GAMMA
            BESTK = K
            BESTWEIGHT = WEIGHT
            BESTPERM = PERM
            print('permutation: {}'.format(PERM))
            print('K={}'.format(K))
            print(PERMSTD)
            print('MX.shape = {}'.format(MX.shape))
            print('PUNCTMATX.shape = {}'.format(PUNCTMATX.shape))
            print('PUNCTSTABX.shape = {}'.format(PUNCTSTABX.shape))
            print('PUNCTLOGX.shape = {}'.format(PUNCTLOGX.shape))
            print('PUNCTSTABZ.shape = {}'.format(PUNCTSTABZ.shape))
            # print('KERZ.shape = {}'.format(KERZ.shape))
            # print('KERX.shape = {}'.format(KERX.shape))
            print('REALLOGX len = {}'.format(KP))
            print('PUNCTLOGZ len = {}'.format(len(PUNCTLOGZ)))
            print('Low weight logical of weight: {}'.format(WEIGHT))
            print('gamma = {}'.format(GAMMA), flush=True)
        (is_dist, _) = distance_lower_bound(PUNCTSTABX, PUNCTLOGX, 3)
        print('Is distance >3: {}'.format(is_dist))
