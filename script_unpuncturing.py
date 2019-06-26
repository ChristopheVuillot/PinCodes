#!/usr/bin/env python3
"""script for unpuncturing codes
"""
# import gc
import numpy as np
import flinalg as fl
from QuantumCodeAnalysis.QuantumCodeAnalysis import low_weight_logical, logicals, get_partner, gram_schmidt, list_low_log
from QuantumCodeAnalysis.puncturing import puncture
from QuantumCodeAnalysis.unpuncturing import unpuncture
from QuantumCodeConstruction.utils import readsparsematrix, writesparsematrix
# from QuantumCodeConstruction.hypergraphproduct import hypergraphproduct
from permutations import PERM2222244, K2222244


def select_unlogx(lowzs, logx):
    anticomm = np.mod(np.dot(lowzs, logx.T), 2)
    print(anticomm)
    nz, nx = anticomm.shape
    numanticom = anticomm.sum(axis=(0))
    notcomplete = numanticom.sum() > 0
    indicesx = []
    while notcomplete:
        indexx = np.argmax(numanticom)
        indicesx.append(indexx)
        for j in range(nz):
            if anticomm[j, indexx] == 1:
                anticomm[j, :] = np.zeros((nx), dtype='uint8')
        numanticom = anticomm.sum(axis=(0))
        notcomplete = numanticom.sum() > 0
    unlogx = []
    leftlogx = []
    for j in indicesx:
        unlogx.append(logx[j, :])
    for j in [k for k in range(nx) if k not in indicesx]:
        leftlogx.append(logx[j, :])
    return np.vstack(unlogx), np.vstack(leftlogx)


PERM = PERM2222244
KINI = K2222244

# MX = readsparsematrix('PCMatrices/systematichp/systematic33_dim3_transpose_7_X.sms').todense()
# MZ = readsparsematrix('PCMatrices/systematichp/systematic33_dim3_transpose_7_Z.sms').todense()
MX = readsparsematrix('PCMatrices/narrowCC/narrowCC_2222244_dim6_X.sms').todense()
MZ = readsparsematrix('PCMatrices/narrowCC/narrowCC_2222244_dim6_Z.sms').todense()
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
print(np.mod(np.dot(PUNCTLOGX, PUNCTLOGX.T), 8))
PUNCTSTABZ = fl.kernel(PUNCTMATX.transpose())
PUNCTK, PUNCTNQ = PUNCTLOGX.shape

# writesparsematrix(PUNCTSTABX, 'PCMatrices/punctured/{}_{}_{}_{}_X.sms'.format('3322233', PUNCTNQ, PUNCTK, 4))
# writesparsematrix(PUNCTSTABZ, 'PCMatrices/punctured/{}_{}_{}_{}_Z.sms'.format('3322233', PUNCTNQ, PUNCTK, 4))

CODES = [(PUNCTSTABX, PUNCTSTABZ, PUNCTLOGX, 4, 'extensivetest_2222244')]


TRIALS = 50
CONTINUE = True

for MX, MZ, LOGXINIT, DINIT, FILENAME in CODES:
    while CONTINUE:
        RX, NQ = MX.shape
        RZ, _ = MZ.shape

        KERX = fl.kernel(MX.T)
        # LOGX, LOGZ = logicals(MX, MZ)
        # LOGX = np.vstack(LOGX)
        LOGX = np.copy(LOGXINIT)

        K = LOGX.shape[0]

        print(NQ, LOGX.shape[0], RZ, RX, K)

        # LOWZS = []

        LOWWEIGHTLOGZ, perm = low_weight_logical(KERX, LOGX, TRIALS)
        # invperm = fl.invert_permutation(perm)
        # LOWWEIGHTLOGZ = LOWWEIGHTLOGZ[invperm]
        # LOWZS.append(LOWWEIGHTLOGZ)
        D = np.count_nonzero(LOWWEIGHTLOGZ)
        print(D)
        # PARTNERIND = get_partner(LOWWEIGHTLOGZ, LOGX)[0]
        # PARTNER = LOGX[PARTNERIND, :]
        # print(PARTNER)
        # LOGX = gram_schmidt(LOWWEIGHTLOGZ, PARTNER, LOGX)
        # #  LOGX = np.delete(LOGX,(PARTNERIND),axis=0)
        # print(PARTNER)
        # UNLOGX = [PARTNER]

        # print('n = {}, k = {}'.format(NQ, LOGX.shape[0]))
        # print('D: {}'.format(D))
        # i = 0
        # while True:
        #     i += 1
        #     print('found {}'.format(i))
        #     LOWWEIGHTLOGZ, perm = low_weight_logical(KERX, LOGX, TRIALS)
        #     invperm = fl.invert_permutation(perm)
        #     LOWWEIGHTLOGZ = LOWWEIGHTLOGZ[invperm]
        #     LOWZS.append(LOWWEIGHTLOGZ)
        #     d = np.count_nonzero(LOWWEIGHTLOGZ)
        #     print(d)
        #     if d > D:
        #         break
        #     PARTNERS = get_partner(LOWWEIGHTLOGZ, LOGX)
        #     if len(PARTNERS):
        #         PARTNERIND = PARTNERS[0]
        #         PARTNER = LOGX[PARTNERIND, :][:]
        #         LOGX = gram_schmidt(LOWWEIGHTLOGZ, PARTNER, LOGX)
        # #      LOGX = np.delete(LOGX,(PARTNERIND),axis=0)
        #         UNLOGX += [PARTNER]
        #         print(d)
        #         print(PARTNERIND)
        # #     print(UNLOGX)
        #     else:
        #         break

        LOWZS = list_low_log(MX, LOGXINIT, DINIT)
        LOWZS = np.vstack(LOWZS)
        UNLOGX, LEFTLOGX = select_unlogx(LOWZS, LOGXINIT)

        LEFTK, _ = LEFTLOGX.shape

        print(UNLOGX.shape)
        LEFTX = np.block([np.zeros((LEFTK, UNLOGX.shape[0]), dtype='uint8'), LEFTLOGX])
        print(LEFTX.shape)

        UNPUNCTMATX = unpuncture(MX, UNLOGX)
        UNPUNCTMATZ = fl.kernel(np.hstack([UNPUNCTMATX.T, LEFTX.T]))

        _, UNPUNCTNQ = UNPUNCTMATX.shape

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

        print('NQ = {}, K = {}, D = {},  nq = {}, k = {}, dx = {}, dz = {}'.format(NQ, K, D, UNPUNCTNQ, UNPUNCTK, UNPUNCTDX, UNPUNCTDZ))
        print('gamma = {}'.format(np.log(UNPUNCTNQ / UNPUNCTK) / np.log(UNPUNCTDZ)))
        if UNPUNCTDZ > D:
            CONTINUE = False
        else:
            MX = UNPUNCTMATX
            MZ = UNPUNCTMATZ
            LOGXINIT = LEFTX
            print(np.mod(np.dot(LEFTX, LEFTX.T), 8))
    writesparsematrix(UNPUNCTMATX, 'PCMatrices/punctured/{}_{}_{}_{}_X.sms'.format(FILENAME, UNPUNCTNQ, UNPUNCTK, UNPUNCTDZ))
    writesparsematrix(UNPUNCTMATZ, 'PCMatrices/punctured/{}_{}_{}_{}_Z.sms'.format(FILENAME, UNPUNCTNQ, UNPUNCTK, UNPUNCTDZ))
    writesparsematrix(LEFTX, 'PCMatrices/punctured/{}_{}_{}_{}_LOGX.sms'.format(FILENAME, UNPUNCTNQ, UNPUNCTK, UNPUNCTDZ))
