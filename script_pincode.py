"""testing script for chaincode.py
"""
# import time
# import scipy.sparse as sp
# import numpy as np
import QuantumCodeConstruction.pincode as pinco
# import QuantumCodeConstruction.hypergraphproduct as hp
# from QuantumCodeConstruction.utils import writesparsematrix, readsparsematrix


if __name__ == "__main__":
    # TRANSITIONS = hp.randomhypergraphproductlist(4, 5, 2, 2, seed=None)
    TRANSTEST = [[[0, 1, 1, 1],
                  [1, 1, 1, 0],
                  [1, 0, 0, 1]],
                 [[1, 0, 1],
                  [0, 1, 1],
                  [1, 1, 0],
                  [1, 0, 1]]]
    POSETTEST = pinco.GrPoset(TRANSTEST, iscomplete=False)
    PTX, PTZ = pinco.pincode(POSETTEST, 1, 1)
    print(POSETTEST.get_flags())
    print('\n')
    PINNEDSETS1 = POSETTEST.get_all_pinned_sets(1)
    for pinset in PINNEDSETS1:
        print(pinset)
        print('\n')
    print(POSETTEST.levelsizes)
    PTCHECKSX, PTQUBITS = PTX.shape
    PTCHECKSZ, _ = PTZ.shape
    print('qubits - xchecks - zchecks = logicals: \n \
           {} - {} - {} = {}'.format(PTQUBITS,
                                     PTCHECKSX,
                                     PTCHECKSZ,
                                     PTQUBITS - PTCHECKSX - PTCHECKSZ))
    # print(PTX)
    # print(PTZ)
    # TRANSITIONS = [readsparsematrix('BoundaryMaps/535_6840_'+str(j)+'.npz').todense()
    #                for j in range(3, 0, -1)]
    # POSETHP = pinco.GrPoset(TRANSITIONS, iscomplete=False)
    # PCX, PCZ = pinco.pincode(POSETHP, 1, 2)
    # print(POSETHP.levelsizes)
    # PCCHECKSX, PCQUBITS = PCX.shape
    # PCCHECKSZ, _ = PCZ.shape
    # print('qubits - xchecks - zchecks = logicals: \n \
    #        {} - {} - {} = {}'.format(PCQUBITS,
    #                                  PCCHECKSX,
    #                                  PCCHECKSZ,
    #                                  PCQUBITS - PCCHECKSX - PCCHECKSZ))
    # # TIME = time.localtime()
    # writesparsematrix(PCX, 'PCMatrices/535_6840_X.sms')
    # writesparsematrix(PCZ, 'PCMatrices/535_6840_Z.sms')
