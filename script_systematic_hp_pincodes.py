"""creating pin codes from systematic hypergraph products
"""
import numpy as np
import QuantumCodeConstruction.pincode as pinco
import QuantumCodeConstruction.hypergraphproduct as hp
from QuantumCodeConstruction.utils import writesparsematrix


if __name__ == "__main__":
    BINARYMATRICESLIST = hp.systematicclassco(4, 3)
    N = len(BINARYMATRICESLIST)
    for j in range(N):
        # the real index of the matrix is j+1 as we skip 0 in the matrix list
        print('Code number {}'.format(j + 1))
        TRANSITIONS = hp.reapeatedhypergraphproduct(BINARYMATRICESLIST[j], 2, transpose=True)
        POSETHP = pinco.GrPoset(TRANSITIONS, iscomplete=False)
        PCX, PCZ = pinco.pincode(POSETHP, 1, 2)
        for k, bmap in enumerate(TRANSITIONS):
            writesparsematrix(bmap, 'BoundaryMaps/systematichp/systematic43_repeat3_transpose_{}_{}.sms'.format(j + 1, k))
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

        CSSCOND = (np.dot(PCX.todense(), PCZ.transpose().todense()) % 2).sum() == 0
        print('This is a valid CSS code: {}'.format(CSSCOND))

        writesparsematrix(PCX, 'PCMatrices/systematichp/systematic43_repeat3_transpose_{}_X.sms'.format(j + 1))
        writesparsematrix(PCZ, 'PCMatrices/systematichp/systematic43_repeat3_transpose_{}_Z.sms'.format(j + 1))
