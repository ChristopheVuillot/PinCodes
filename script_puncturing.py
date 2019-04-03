"""script for puncturing codes
"""
# import gc
import numpy as np
import flinalg as fl
from QuantumCodeAnalysis.QuantumCodeAnalysis import low_weight_logical, distance_lower_bound
from QuantumCodeAnalysis.puncturing import puncture
from QuantumCodeConstruction.utils import readsparsematrix


# MX = readsparsematrix('PCMatrices/systematichp/systematic33_dim3_transpose_10_X.sms').todense()
# MZ = readsparsematrix('PCMatrices/systematichp/systematic33_dim3_transpose_10_Z.sms').todense()
# MX = readsparsematrix('PCMatrices/narrowCC/narrowCC2_dim6_X.sms').todense()
# MZ = readsparsematrix('PCMatrices/narrowCC/narrowCC2_dim6_Z.sms').todense()
MX = readsparsematrix('PCMatrices/narrowCC/narrowCC_2222233_dim6_X.sms').todense()
MZ = readsparsematrix('PCMatrices/narrowCC/narrowCC_2222233_dim6_Z.sms').todense()
# MX = readsparsematrix('PCMatrices/535_3420_XCOS.sms').todense()
# MZ = readsparsematrix('PCMatrices/535_3420_ZCOS.sms').todense()


RX, NQ = MX.shape
RZ, _ = MZ.shape

# MX = np.array(MX, dtype='uint8')
# SX = fl.row_reduce_transform(MX)

BESTGAMMA = 5
for k in range(17, 18):
    print('trying {} punctures:'.format(k), flush=True)
    for _ in range(1):
        # PERM = np.random.permutation(NQ)
        PERM = [ 55,  75, 103,  41, 172, 135, 162, 156,  24, 101, 113,  50,  73, 120, 148,  58,   9,  85,
                 81,  12,  51, 182,  67,  14, 123, 178, 174, 122,  34, 104,  39, 164, 191,  70,  54,   8,
                 10, 126, 151, 131,   3,  42,  56, 149, 102,  65,  31,  57, 108,  62, 109, 134, 165,  47,
                 52,  40,  17, 187, 144,  79, 154,  90,  78,  26,  30,  19, 176,   4, 107,  49, 140,  28,
                  2, 145, 163,  45, 106, 179, 115, 175, 110,  44,  53,  15, 128,  93, 158,  33,  84,  37,
                 60,  46, 180,  23,  76, 111,  35,  87, 142,  83,  95, 169, 185,  92, 166,  32, 147, 160,
                137,   5,  22, 136, 177, 121, 105, 150,  88, 161, 138, 170,  29,  38, 146, 186, 171,  25,
                189, 116,  11, 183, 114, 188,  82,  61,  94,  64,  66, 141, 157,  72, 133,  89,  43, 125,
                 16, 124, 129,  86,   7,  27,  69, 130,  80,  71, 100, 173, 119,  68, 184,  20,  18,  97,
                 91,  13, 117, 152,  48,  36, 143, 181, 167, 132,   1,  96,  99, 127,   6, 139, 112,   0,
                 63,  59, 190, 118,  74, 168, 155,  77, 153,  98,  21, 159]
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
        (is_dist_4, _) = distance_lower_bound(PUNCTSTABX, PUNCTLOGX, 3)
        print('Is distance 4: {}'.format(is_dist_4))
