"""creating pin codes from narrow chain complex
"""
import numpy as np
import QuantumCodeConstruction.pincode as pinco
from QuantumCodeConstruction.utils import writesparsematrix


if __name__ == "__main__":
    D = 6
    X = int(D/3)
    Z = D-X
    SIZE = 2
    # print('Generating narrow Pin Code with D={}, x={}, z={}, q={}'.format(D, X, Z, SIZE))
    # TRANSITIONS = [np.ones((SIZE, SIZE), dtype='uint8') for _ in range(D)]
    print('Narrow pin code 2-2-2-2-2-4-sep-4')
    M44 = np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 1, 1], [0, 0, 1, 1]], dtype='uint8')
    TRANSITIONS = [np.ones((2, 2), dtype='uint8'),
                   np.ones((2, 2), dtype='uint8'),
                   np.ones((2, 2), dtype='uint8'),
                   np.ones((2, 2), dtype='uint8'),
                   np.ones((2, 4), dtype='uint8'),
                   M44]
    POSETHP = pinco.GrPoset(TRANSITIONS, iscomplete=False)
    PCX, PCZ = pinco.pincode(POSETHP, X, Z)
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

    writesparsematrix(PCX, 'PCMatrices/narrowCC/narrowCC_222224sep4_dim{}_X.sms'.format(D))
    writesparsematrix(PCZ, 'PCMatrices/narrowCC/narrowCC_222224sep4_dim{}_Z.sms'.format(D))
