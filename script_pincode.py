"""testing script for chaincode.py
"""
# import time
# import scipy.sparse as sp
# import numpy as np
from itertools import combinations
from scipy.sparse import dok_matrix
import QuantumCodeConstruction.pincode as pinco
# import QuantumCodeConstruction.hypergraphproduct as hp
from QuantumCodeConstruction.utils import readsparsematrix, writesparsematrix


def readcosetfile(filename):
    """read file with one coset per line
    """
    with open(filename, 'r') as cosetfile:
        return [[int(i) for i in co.split(',')] for co in cosetfile.readlines()]

if __name__ == "__main__":
    # TRANSITIONS = hp.randomhypergraphproductlist(4, 5, 2, 2, seed=None)
    #  TRANSTEST = [[[0, 1, 1, 1],
    #                [1, 1, 1, 0]],
    #               [[1, 0],
    #                [0, 1],
    #                [1, 1],
    #                [1, 0]]]
    # TRANSTEST = [np.ones([3, 4]), np.ones([4, 4]), np.ones([4, 3])]
    # POSETTEST = pinco.GrPoset(TRANSTEST, iscomplete=False)
    # PTX, PTZ = pinco.pincode(POSETTEST, 1, 1)
    # print(POSETTEST.get_flags())
    # print('\n')
    # PINNEDSETS1 = POSETTEST.get_all_pinned_sets(1)
    # PINNEDSETS2 = POSETTEST.get_all_pinned_sets(2)
    # for pinset in PINNEDSETS1:
    #     print(pinset)
    #     print('\n')
    # print(len(PINNEDSETS1))
    # print(len(PINNEDSETS2))
    # print({len(a) for a in PINNEDSETS1})
    # print({len(a) for a in PINNEDSETS2})
    # print(POSETTEST.levelsizes)
    # PTCHECKSX, PTQUBITS = PTX.shape
    # PTCHECKSZ, _ = PTZ.shape
    # print('qubits - xchecks - zchecks = logicals: \n \
    #        {} - {} - {} = {}'.format(PTQUBITS,
    #                                  PTCHECKSX,
    #                                  PTCHECKSZ,
    #                                  PTQUBITS - PTCHECKSX - PTCHECKSZ))
    TRANSITIONS = [readsparsematrix('BoundaryMaps/535_6840_'+str(j)+'.npz').todense()
                   for j in range(3, 0, -1)]
    POSETHP = pinco.GrPoset(TRANSITIONS, iscomplete=False)
    PCX, PCZ = pinco.pincode(POSETHP, 1, 2)
    WEIGHTSX = {a for che in PCX.sum(axis=1).tolist() for a in che}
    WEIGHTSZ = {a for che in PCZ.sum(axis=1).tolist() for a in che}
    print(WEIGHTSX)
    print(WEIGHTSZ)
    print(POSETHP.levelsizes)
    PCCHECKSX, PCQUBITS = PCX.shape
    PCCHECKSZ, _ = PCZ.shape
    print('qubits - xchecks - zchecks = logicals: \n \
           {} - {} - {} = {}'.format(PCQUBITS,
                                     PCCHECKSX,
                                     PCCHECKSZ,
                                     PCQUBITS - PCCHECKSX - PCCHECKSZ))
    ALL1COSETS = []
    for j in range(4):
        ALL1COSETS += readcosetfile('PCMatrices/535_6840_'+str(j)+'.cos')
    ALL2COSETS = []
    for cos1, cos2 in combinations(ALL1COSETS, 2):
        inter = set(cos1).intersection(set(cos2))
        if inter:
            ALL2COSETS.append(list(inter))
    print(len(ALL1COSETS))
    print({len(co) for co in ALL1COSETS})
    print(len(ALL2COSETS))
    print({len(co) for co in ALL2COSETS})
    NQ = max([a for sub in ALL1COSETS for a in sub])
    print(NQ)
    PCXCOS = dok_matrix((len(ALL1COSETS), NQ), dtype='int')
    for j, xche in enumerate(ALL1COSETS):
        for elem in xche:
            PCXCOS[j, elem-1] = 1
    PCZCOS = dok_matrix((len(ALL2COSETS), NQ), dtype='int')
    for j, zche in enumerate(ALL2COSETS):
        for elem in zche:
            PCZCOS[j, elem-1] = 1
    # TIME = time.localtime()
    writesparsematrix(PCX, 'PCMatrices/535_6840_X.sms')
    writesparsematrix(PCZ, 'PCMatrices/535_6840_Z.sms')
    writesparsematrix(PCXCOS, 'PCMatrices/535_6840_XCOS.sms')
    writesparsematrix(PCZCOS, 'PCMatrices/535_6840_ZCOS.sms')
