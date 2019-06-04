#!/usr/bin/env python3
"""script for unpuncturing codes
"""
# import gc
import numpy as np
import flinalg as fl
from QuantumCodeAnalysis.QuantumCodeAnalysis import low_weight_logical, distance_lower_bound, logicals, get_partner, gram_schmidt
from QuantumCodeAnalysis.unpuncturing import unpuncture
from QuantumCodeConstruction.utils import readsparsematrix
from QuantumCodeConstruction.hypergraphproduct import hypergraphproduct
# from permutations import PERM2222244, K2222244


MX = readsparsematrix('PCMatrices/systematichp/systematic33_dim3_transpose_7_X.sms').todense()
MZ = readsparsematrix('PCMatrices/systematichp/systematic33_dim3_transpose_7_Z.sms').todense()
# MX = readsparsematrix('PCMatrices/narrowCC/narrowCC2_dim6_X.sms').todense()
# MZ = readsparsematrix('PCMatrices/narrowCC/narrowCC2_dim6_Z.sms').todense()
# MX = readsparsematrix('PCMatrices/narrowCC/narrowCC_2222244_dim6_X.sms').todense()
# MZ = readsparsematrix('PCMatrices/narrowCC/narrowCC_2222244_dim6_Z.sms').todense()
#MX = readsparsematrix('PCMatrices/535_3420_XCOS.sms').todense()
#MZ = readsparsematrix('PCMatrices/535_3420_ZCOS.sms').todense()

a = np.array(np.vstack([np.roll(np.array([1,1,0,0],dtype='uint8'),i) for i in range(4)]),dtype='uint8')
b = np.array(np.vstack([np.roll(np.array([1,1,0],dtype='uint8'),i) for i in range(3)]),dtype='uint8')
MX,MZ = hypergraphproduct(a,b)
MX = np.mod(np.array(MX,dtype='uint8'),2)
MZ = np.mod(np.array(MZ,dtype='uint8'),2)

CODES = [(MX,MZ)]


TRIALS = 30

for MX,MZ in CODES:
    RX, NQ = MX.shape
    RZ, _ = MZ.shape

    KERZ = fl.kernel(MZ.T)
    LOGX, LOGZ = logicals(MX,MZ)
    LOGX = np.vstack(LOGX)

    K = LOGX.shape[0]

    print(NQ,LOGX.shape[0],RZ,RX)

    LOWWEIGHTLOGZ, perm = low_weight_logical(KERZ, LOGX, TRIALS)
    invperm = fl.invert_permutation(perm)
    LOWWEIGHTLOGZ = LOWWEIGHTLOGZ[invperm]
    D = np.count_nonzero(LOWWEIGHTLOGZ)
    print(D)
    PARTNERIND = get_partner(LOWWEIGHTLOGZ,LOGX)[0]
    PARTNER = LOGX[PARTNERIND,:]
    print(PARTNER)
    LOGX = gram_schmidt(PARTNER,LOGX)
#    LOGX = np.delete(LOGX,(PARTNERIND),axis=0)
    print(PARTNER)
    UNLOGX = [PARTNER]

    print('n = {}, k = {}'.format(NQ,LOGX.shape[0]))
    print('D: {}'.format(D))
    i=0
    while True:
        i += 1
        print('found {}'.format(i))
        LOWWEIGHTLOGZ, perm = low_weight_logical(KERZ, LOGX, TRIALS)
        invperm = fl.invert_permutation(perm)
        LOWWEIGHTLOGZ = LOWWEIGHTLOGZ[invperm]
        d = np.count_nonzero(LOWWEIGHTLOGZ)
        print(d)
        if d > D:
            break
        PARTNERS = get_partner(LOWWEIGHTLOGZ,LOGX)
        if PARTNERS:
            PARTNERIND = PARTNERS[0]
            PARTNER = LOGX[PARTNERIND,:][:]
            LOGX = gram_schmidt(PARTNER,LOGX)
#        LOGX = np.delete(LOGX,(PARTNERIND),axis=0)
            UNLOGX += [PARTNER]
            print(d)
            print(PARTNERIND)
#        print(UNLOGX)
        else:
            break

    UNLOGX = np.vstack(UNLOGX)
    LEFTX = np.block([np.zeros((K,UNLOGX.shape[0]),dtype='uint8'),LOGX])

    UNPUNCTMATX = unpuncture(MX, UNLOGX)
    UNPUNCTMATZ = fl.kernel(np.hstack([UNPUNCTMATX.T,LEFTX.T]))

    print(MX)
    print(MZ)
    print(UNPUNCTMATX)
    print(UNPUNCTMATZ)

    UNPUNCTKERZ = fl.kernel(UNPUNCTMATZ.T)
    UNPUNCTLOGX = fl.quotient_basis(UNPUNCTKERZ,UNPUNCTMATX)
    UNPUNCTLOGX = np.vstack(UNPUNCTLOGX)

    UNPUNCTLOWWEIGHTLOGZ, _ = low_weight_logical(UNPUNCTKERZ, UNPUNCTLOGX, TRIALS)
    print(UNPUNCTLOWWEIGHTLOGZ)
    UNPUNCTD = np.count_nonzero(UNPUNCTLOWWEIGHTLOGZ)

    print('D = {},  d = {}'.format(D,UNPUNCTD))
