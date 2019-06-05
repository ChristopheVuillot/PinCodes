#!/usr/bin/env python3
"""script for unpuncturing codes
"""
# import gc
import numpy as np
import flinalg as fl
from QuantumCodeAnalysis.QuantumCodeAnalysis import low_weight_logical, logicals, get_partner, gram_schmidt
from QuantumCodeAnalysis.puncturing import puncture
from QuantumCodeAnalysis.unpuncturing import unpuncture
from QuantumCodeConstruction.utils import readsparsematrix
# from QuantumCodeConstruction.hypergraphproduct import hypergraphproduct
from permutations import PERM3322233, K3322233

PERM = PERM3322233
KINI = K3322233

# MX = readsparsematrix('PCMatrices/systematichp/systematic33_dim3_transpose_7_X.sms').todense()
# MZ = readsparsematrix('PCMatrices/systematichp/systematic33_dim3_transpose_7_Z.sms').todense()
MX = readsparsematrix('PCMatrices/narrowCC/narrowCC_3322233_dim6_X.sms').todense()
MZ = readsparsematrix('PCMatrices/narrowCC/narrowCC_3322233_dim6_Z.sms').todense()
# MX = readsparsematrix('PCMatrices/narrowCC/narrowCC_2222244_dim6_X.sms').todense()
# MZ = readsparsematrix('PCMatrices/narrowCC/narrowCC_2222244_dim6_Z.sms').todense()
# MX = readsparsematrix('PCMatrices/535_3420_XCOS.sms').todense()
# MZ = readsparsematrix('PCMatrices/535_3420_ZCOS.sms').todense()

# a = np.array(np.vstack([np.roll(np.array([1,1,0,0],dtype='uint8'),i) for i in range(4)]),dtype='uint8')
# b = np.array(np.vstack([np.roll(np.array([1,1,0],dtype='uint8'),i) for i in range(3)]),dtype='uint8')
# MX,MZ = hypergraphproduct(a,b)
# MX = np.mod(np.array(MX,dtype='uint8'),2)
# MZ = np.mod(np.array(MZ,dtype='uint8'),2)

PUNCTMATX, PERMSTD = puncture(MX, PERM, KINI)
PUNCTSTABX = np.array(PUNCTMATX[KINI:, :], dtype='uint8')
PUNCTLOGX = np.array(PUNCTMATX[:KINI, :], dtype='uint8')
PUNCTSTABZ = fl.kernel(PUNCTMATX.transpose())


CODES = [(PUNCTSTABX, PUNCTSTABZ)]


TRIALS = 50

for MX, MZ in CODES:
    RX, NQ = MX.shape
    RZ, _ = MZ.shape

    KERX = fl.kernel(MX.T)
    LOGX, LOGZ = logicals(MX, MZ)
    LOGX = np.vstack(LOGX)

    K = LOGX.shape[0]

    print(NQ, LOGX.shape[0], RZ, RX, K)

    LOWWEIGHTLOGZ, perm = low_weight_logical(KERX, LOGX, TRIALS)
    invperm = fl.invert_permutation(perm)
    LOWWEIGHTLOGZ = LOWWEIGHTLOGZ[invperm]
    D = np.count_nonzero(LOWWEIGHTLOGZ)
    print(D)
    PARTNERIND = get_partner(LOWWEIGHTLOGZ, LOGX)[0]
    PARTNER = LOGX[PARTNERIND, :]
    print(PARTNER)
    LOGX = gram_schmidt(LOWWEIGHTLOGZ, PARTNER, LOGX)
#    LOGX = np.delete(LOGX,(PARTNERIND),axis=0)
    print(PARTNER)
    UNLOGX = [PARTNER]

    print('n = {}, k = {}'.format(NQ, LOGX.shape[0]))
    print('D: {}'.format(D))
    i = 0
    while True:
        i += 1
        print('found {}'.format(i))
        LOWWEIGHTLOGZ, perm = low_weight_logical(KERX, LOGX, TRIALS)
        invperm = fl.invert_permutation(perm)
        LOWWEIGHTLOGZ = LOWWEIGHTLOGZ[invperm]
        d = np.count_nonzero(LOWWEIGHTLOGZ)
        print(d)
        if d > D:
            break
        PARTNERS = get_partner(LOWWEIGHTLOGZ, LOGX)
        if len(PARTNERS):
            PARTNERIND = PARTNERS[0]
            PARTNER = LOGX[PARTNERIND, :][:]
            LOGX = gram_schmidt(LOWWEIGHTLOGZ, PARTNER, LOGX)
#        LOGX = np.delete(LOGX,(PARTNERIND),axis=0)
            UNLOGX += [PARTNER]
            print(d)
            print(PARTNERIND)
#        print(UNLOGX)
        else:
            break

    UNLOGX = np.vstack(UNLOGX)
    LEFTX = np.block([np.zeros((K, UNLOGX.shape[0]), dtype='uint8'), LOGX])

    UNPUNCTMATX = unpuncture(MX, UNLOGX)
    UNPUNCTMATZ = fl.kernel(np.hstack([UNPUNCTMATX.T, LEFTX.T]))

    CSSCOND = np.mod(np.dot(UNPUNCTMATX, UNPUNCTMATZ.T), 2).sum() == 0
    print(CSSCOND)

    # print(MX)
    # print(MZ)
    # print(UNPUNCTMATX)
    # print(UNPUNCTMATZ)

    UNPUNCTKERX = fl.kernel(UNPUNCTMATX.T)
    UNPUNCTKERZ = fl.kernel(UNPUNCTMATZ.T)
    UNPUNCTLOGX, UNPUNCTLOGZ = logicals(UNPUNCTMATX, UNPUNCTMATZ)
    UNPUNCTK = len(UNPUNCTLOGX)
    UNPUNCTLOGX = np.vstack(UNPUNCTLOGX)

    UNPUNCTLOWWEIGHTLOGZ, _ = low_weight_logical(UNPUNCTKERX, UNPUNCTLOGX, TRIALS)
    UNPUNCTLOWWEIGHTLOGX, _ = low_weight_logical(UNPUNCTKERZ, UNPUNCTLOGZ, TRIALS)
    # print(UNPUNCTLOWWEIGHTLOGZ)
    UNPUNCTDZ = np.count_nonzero(UNPUNCTLOWWEIGHTLOGZ)
    UNPUNCTDX = np.count_nonzero(UNPUNCTLOWWEIGHTLOGX)

    print('D = {},  dx = {}, dz = {}, K = {}, k = {}'.format(D, UNPUNCTDX, UNPUNCTDZ, K, UNPUNCTK))
