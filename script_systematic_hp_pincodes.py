"""creating pin codes from systematic hypergraph products
"""
import numpy as np
import QuantumCodeConstruction.pincode as pinco
import QuantumCodeConstruction.hypergraphproduct as hp
from QuantumCodeConstruction.utils import writesparsematrix


if __name__ == "__main__":
    D = 2
    X = 1  # int(D / 3)
    Z = D - X
    CHECKS = 3
    BITS = 4
    BINARYMATRICESLIST = hp.systematicclassco(CHECKS, BITS)
    N = len(BINARYMATRICESLIST)
    print('Generating Pin Codes with D={}, x={}, z={}, from all {}x{} binary matrices'.format(D, X, Z, CHECKS, BITS))
    for j in range(N):
        # the real index of the matrix is j+1 as we skip 0 in the matrix list
        print('Code systematic{}{}_dim{}_transpose_{}'.format(CHECKS, BITS, D, j + 1))
        TRANSITIONS = hp.reapeatedhypergraphproduct(BINARYMATRICESLIST[j], D - 1, transpose=True)
        POSETHP = pinco.GrPoset(TRANSITIONS, iscomplete=False)
        PCX, PCZ = pinco.pincode(POSETHP, X, Z)
        for k, bmap in enumerate(TRANSITIONS):
            writesparsematrix(bmap, 'BoundaryMaps/systematichp/systematic{}{}_dim{}_transpose_{}_{}.sms'.format(CHECKS, BITS, D, j + 1, k))
        # WEIGHTSX = {a for che in PCX.sum(axis=1).tolist() for a in che}
        # WEIGHTSZ = {a for che in PCZ.sum(axis=1).tolist() for a in che}
        # print(WEIGHTSX)
        # print(WEIGHTSZ)
        print(POSETHP.levelsizes)
        PCCHECKSX, PCQUBITS = PCX.shape
        PCCHECKSZ, _ = PCZ.shape
        print('qubits - xchecks - zchecks = logicals: \n \
               {} - {} - {} = {}'.format(PCQUBITS,
                                         PCCHECKSX,
                                         PCCHECKSZ,
                                         PCQUBITS - PCCHECKSX - PCCHECKSZ))

        CSSCOND = (np.dot(PCX.todense(), PCZ.transpose().todense()) % 2).sum() == 0
        print('This is a valid CSS code: {}'.format(CSSCOND))

        writesparsematrix(PCX, 'PCMatrices/systematichp/systematic{}{}_dim{}_transpose_{}_X.sms'.format(CHECKS, BITS, D, j + 1))
        writesparsematrix(PCZ, 'PCMatrices/systematichp/systematic{}{}_dim{}_transpose_{}_Z.sms'.format(CHECKS, BITS, D, j + 1))
