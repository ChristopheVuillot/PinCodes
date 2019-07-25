"""creating CC..CZ pin codes from narrow chain complex
using N. Rengaswamy technique (for Reed-Muller codes)
adapted to pin codes
"""
import numpy as np
import QuantumCodeConstruction.pincode as pinco
from QuantumCodeConstruction.utils import writesparsematrix


if __name__ == "__main__":
    D = 5
    LX = 2
    X = LX - 1
    print('Generating CC..CZ Pin Code 222224')
    M44 = np.array([[1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 1, 1], [1, 0, 0, 1]], dtype='uint8')
    M33 = np.array([[1, 1, 0], [0, 1, 1], [1, 0, 1]], dtype='uint8')
    TRANSITIONS = [np.ones((2, 2), dtype='uint8'),
                   np.ones((2, 2), dtype='uint8'),
                   np.ones((2, 2), dtype='uint8'),
                   np.ones((2, 2), dtype='uint8'),
                   np.ones((2, 4), dtype='uint8')]
    POSETHP = pinco.GrPoset(TRANSITIONS, iscomplete=False)
    PCX, PCZ = pinco.cczpincode(POSETHP, LX)
    print(PCX)
    print(PCZ)
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

    CSSCOND = (np.dot(PCX, PCZ.transpose()) % 2).sum() == 0
    print('This is a valid CSS code: {}'.format(CSSCOND))

    writesparsematrix(PCX, 'PCMatrices/notgauge/cczpinco_222224_dim{}_X.sms'.format(D))
    writesparsematrix(PCZ, 'PCMatrices/notgauge/cczpinco_222224_dim{}_Z.sms'.format(D))
