"""creating pin codes from systematic hypergraph products
"""
import numpy as np
import QuantumCodeConstruction.pincode as pinco
import QuantumCodeConstruction.hypergraphproduct as hp
from QuantumCodeConstruction.utils import writesparsematrix


if __name__ == "__main__":
    BINARYMATRICESLIST = hp.systematicclassco(3, 3)
    N = len(BINARYMATRICESLIST)
    for j in range(N):
        TRANSITIONS = hp.reapeatedhypergraphproduct(BINARYMATRICESLIST[j], 2, transpose=True)
        POSETHP = pinco.GrPoset(TRANSITIONS, iscomplete=False)
        PCX, PCZ = pinco.pincode(POSETHP, 1, 2)
        for k, bmap in enumerate(TRANSITIONS):
            writesparsematrix(bmap, 'BoundaryMaps/systematichp/systematic33_repeat3_transpose_{}_{}.sms'.format(j, k))
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

        writesparsematrix(PCX, 'PCMatrices/systematichp/systematic33_repeat3_transpose_{}_X.sms'.format(j))
        writesparsematrix(PCZ, 'PCMatrices/systematichp/systematic33_repeat3_transpose_{}_Z.sms'.format(j))
