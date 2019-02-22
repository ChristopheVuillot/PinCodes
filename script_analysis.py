"""Script to analyse codes
"""
import numpy as np
from QuantumCodeAnalysis.QuantumCodeAnalysis import low_weight_logical, logicals
from QuantumCodeConstruction.utils import readsparsematrix

for SEED in range(1000, 1100):
    if SEED == 1097:
        continue
    MX = readsparsematrix('PCMatrices/randomhp_{}_X.sms'.format(SEED)).todense()
    MZ = readsparsematrix('PCMatrices/randomhp_{}_Z.sms'.format(SEED)).todense()

    print('Properties of the code randomhp_{}:'.format(SEED))
    print('X-check matrix is {}x{}'.format(MX.shape[0], MX.shape[1]))
    print('Z-check matrix is {}x{}'.format(MZ.shape[0], MZ.shape[1]))

    XW = {MX[j, :].sum() for j in range(MX.shape[0])}
    ZW = {MZ[j, :].sum() for j in range(MZ.shape[0])}

    print('X-checks have weigths: {}'.format(XW))
    print('Z-checks have weigths: {}'.format(ZW))

    (LX, _), (LZ, _) = logicals(MX, MZ)
    if LX:
        LX = np.vstack(LX)
        LZ = np.vstack(LZ)
        GX = np.block([[LX], [MX]])
        GZ = np.block([[LZ], [MZ]])
    else:
        print('Zero logical operator !')
        continue

    print('There are {} logical qubits'.format(LX.shape[0]))

    LOWX, PERMX = low_weight_logical(GX, LZ, 1)
    LOWZ, PERMZ = low_weight_logical(GZ, LX, 1)

    print('There is a X-logical operator of weight {}'.format(LOWX.sum()))
    print('There is a Z-logical operator of weight {}'.format(LOWZ.sum()))

    print('Check correct x log: {}'.format((np.dot(MZ[:, PERMX], LOWX) % 2).sum() == 0))
    print('Check correct z log: {}'.format((np.dot(MX[:, PERMZ], LOWZ) % 2).sum() == 0))
