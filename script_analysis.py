"""Script to analyse codes
"""
from itertools import combinations
import numpy as np
from QuantumCodeAnalysis.QuantumCodeAnalysis import low_weight_logical, logicals
from QuantumCodeConstruction.utils import readsparsematrix

goodones = []
for SEED in range(1, 2):  # 2**9-1):
    # if SEED == 1097:
    #     continue
    # MX = readsparsematrix('PCMatrices/systematichp/systematic33_dim3_transpose_{}_X.sms'.format(SEED)).todense()
    # MZ = readsparsematrix('PCMatrices/systematichp/systematic33_dim3_transpose_{}_Z.sms'.format(SEED)).todense()
    # MX = readsparsematrix('PCMatrices/randomhp/randomhp_swap_{}_X.sms'.format(SEED)).todense()
    # MZ = readsparsematrix('PCMatrices/randomhp/randomhp_swap_{}_Z.sms'.format(SEED)).todense()
    MZ = readsparsematrix('PCMatrices/color_code_535/6840_abcbadcbabadcbabadcbadcbabdcbabdX.sms').todense()
    MX = readsparsematrix('PCMatrices/color_code_535/6840_abcbadcbabadcbabadcbadcbabdcbabdZ.sms').todense()
    # MZ = readsparsematrix('PCMatrices/535_6840_Z.sms').todense()
    # MX = readsparsematrix('PCMatrices/535_6840_X.sms').todense()

    print('Properties of the code 3420_abcdcbadcbabadcbabadcbdcdcb:{}'.format(' '))
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

    print('logical X')
    LOWX, PERMX = low_weight_logical(GX, LZ, 1)
    print('logical Z')
    LOWZ, PERMZ = low_weight_logical(GZ, LX, 1)

    print('There is a X-logical operator of weight {}'.format(LOWX.sum()))
    print('There is a Z-logical operator of weight {}'.format(LOWZ.sum()))

    print('Check correct x log: {}'.format((np.dot(MZ[:, PERMX], LOWX) % 2).sum() == 0))
    print('Check correct z log: {}'.format((np.dot(MX[:, PERMZ], LOWZ) % 2).sum() == 0))

    print('Checking tri-orthogonal condition |L_j wedge L_k wedge S_l|')
    K, N = LX.shape
    NSX, _ = MX.shape
    TRICOND = True
    wrong = []
    for ilog1, ilog2 in combinations(range(K), 2):
        for istab in range(NSX):
            TEST = np.multiply(LX[ilog1, :],
                               np.multiply(LX[ilog2, :],
                                           MX[istab, :])).sum() % 2 == 0
            TRICOND = TRICOND and TEST
            if not TEST:
                wrong.append(((ilog1, ilog2), istab))
            if wrong:
                break
        if wrong:
            break
    print('Triorthogonal: {}'.format(TRICOND))
    if TRICOND:
        goodones.append(SEED)
    print(len(wrong))
# print(goodones)
    # print(LX[wrong[0][0][0], :])
    # print(LX[wrong[0][0][1], :])
    # print(MX[wrong[0][1], :])
    # a = wrong[0][0][0]
    # b = wrong[0][0][1]
    # c = wrong[0][1]
    # count = 0
    # for j in range(N):
    #     count += LX[a, j]*LX[b, j]*MX[c, j]
    #     print(LX[a, j], LX[b, j], MX[c, j])
    # print(count)
