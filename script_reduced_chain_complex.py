"""script to test reduced chain complexes of pin codes
"""
from itertools import combinations
import numpy as np
from QuantumCodeConstruction.pincode import GrPoset, reduced_chain_complex, pincode, reduced_vec_to_flag_vec
from QuantumCodeConstruction.utils import readsparsematrix
from QuantumCodeAnalysis.QuantumCodeAnalysis import logicals  # , low_weight_logical


TRANSITIONS = [readsparsematrix('BoundaryMaps/systematichp/systematic33_dim3_transpose_7_{}.sms'.format(j)).todense() for j in range(3)]
# TRANSITIONS = [readsparsematrix('BoundaryMaps/535_3420_{}.sms'.format(j)).todense().transpose() for j in range(1, 4)]
POSET = GrPoset(TRANSITIONS, iscomplete=False)
print('correct boundary map: {}'.format(POSET.check_boundary_map()))
PCX, PCZ = pincode(POSET, 1, 2)
RX, NQ = PCX.shape
print('correct boundary map: {}'.format(POSET.check_boundary_map()))
(LX, _), (LZ, _) = logicals(PCX.todense(), PCZ.todense())
print(len(LX))
# print(PCX.todense())
# print(PCZ.todense())

print('Constructing reduced Z logicals')
XMAT0, ZMAT0, REDQ0 = reduced_chain_complex(POSET, (0,))
(LX0, _), (LZ0, _) = logicals(XMAT0, ZMAT0)
print(len(LZ0))
FLAGVECS0 = np.zeros((len(LZ0), NQ), dtype='uint8')
for j in range(len(LZ0)):
    FLAGVECS0[j] = reduced_vec_to_flag_vec(POSET, LZ0[j], REDQ0)
print('are actual logicals if 0: {}'.format((np.dot(PCX.todense(), FLAGVECS0.transpose()) % 2).sum()))

XMAT1, ZMAT1, REDQ1 = reduced_chain_complex(POSET, (1,))
(LX1, _), (LZ1, _) = logicals(XMAT1, ZMAT1)
print(len(LZ1))
FLAGVECS1 = np.zeros((len(LZ1), NQ), dtype='uint8')
for j in range(len(LZ1)):
    FLAGVECS1[j] = reduced_vec_to_flag_vec(POSET, LZ1[j], REDQ1)
print('are actual logicals if 0: {}'.format((np.dot(PCX.todense(), FLAGVECS1.transpose()) % 2).sum()))

XMAT2, ZMAT2, REDQ2 = reduced_chain_complex(POSET, (2,))
(LX2, _), (LZ2, _) = logicals(XMAT2, ZMAT2)
print(len(LZ2))
FLAGVECS2 = np.zeros((len(LZ2), NQ), dtype='uint8')
for j in range(len(LZ2)):
    FLAGVECS2[j] = reduced_vec_to_flag_vec(POSET, LZ2[j], REDQ2)
print('are actual logicals if 0: {}'.format((np.dot(PCX.todense(), FLAGVECS2.transpose()) % 2).sum()))

XMAT3, ZMAT3, REDQ3 = reduced_chain_complex(POSET, (3,))
(LX3, _), (LZ3, _) = logicals(XMAT3, ZMAT3)
print(len(LZ3))
FLAGVECS3 = np.zeros((len(LZ3), NQ), dtype='uint8')
for j in range(len(LZ3)):
    FLAGVECS3[j] = reduced_vec_to_flag_vec(POSET, LZ3[j], REDQ3)
print('are actual logicals if 0: {}'.format((np.dot(PCX.todense(), FLAGVECS3.transpose()) % 2).sum()))

print('Constructing reduced X logicals')
XMAT01, ZMAT01, REDQ01 = reduced_chain_complex(POSET, (0, 1))
(LX01, _), (LZ01, _) = logicals(XMAT01, ZMAT01)
print(len(LZ01))
FLAGVECS01 = np.zeros((len(LX01), NQ), dtype='uint8')
for j in range(len(LX01)):
    FLAGVECS01[j] = reduced_vec_to_flag_vec(POSET, LX01[j], REDQ01)
print('are actual logicals if 0: {}'.format((np.dot(PCX.todense(), FLAGVECS01.transpose()) % 2).sum()))

XMAT02, ZMAT02, REDQ02 = reduced_chain_complex(POSET, (0, 2))
(LX02, _), (LZ02, _) = logicals(XMAT02, ZMAT02)
print(len(LZ02))
FLAGVECS02 = np.zeros((len(LX02), NQ), dtype='uint8')
for j in range(len(LX02)):
    FLAGVECS02[j] = reduced_vec_to_flag_vec(POSET, LX02[j], REDQ02)
print('are actual logicals if 0: {}'.format((np.dot(PCX.todense(), FLAGVECS02.transpose()) % 2).sum()))

XMAT03, ZMAT03, REDQ03 = reduced_chain_complex(POSET, (0, 3))
(LX03, _), (LZ03, _) = logicals(XMAT03, ZMAT03)
print(len(LZ03))
FLAGVECS03 = np.zeros((len(LX03), NQ), dtype='uint8')
for j in range(len(LX03)):
    FLAGVECS03[j] = reduced_vec_to_flag_vec(POSET, LX03[j], REDQ03)
print('are actual logicals if 0: {}'.format((np.dot(PCX.todense(), FLAGVECS03.transpose()) % 2).sum()))

XMAT12, ZMAT12, REDQ12 = reduced_chain_complex(POSET, (1, 2))
(LX12, _), (LZ12, _) = logicals(XMAT12, ZMAT12)
print(len(LZ12))
FLAGVECS12 = np.zeros((len(LX12), NQ), dtype='uint8')
for j in range(len(LX12)):
    FLAGVECS12[j] = reduced_vec_to_flag_vec(POSET, LX12[j], REDQ12)
print('are actual logicals if 0: {}'.format((np.dot(PCX.todense(), FLAGVECS12.transpose()) % 2).sum()))

XMAT13, ZMAT13, REDQ13 = reduced_chain_complex(POSET, (1, 3))
(LX13, _), (LZ13, _) = logicals(XMAT13, ZMAT13)
print(len(LZ13))
FLAGVECS13 = np.zeros((len(LX13), NQ), dtype='uint8')
for j in range(len(LX13)):
    FLAGVECS13[j] = reduced_vec_to_flag_vec(POSET, LX13[j], REDQ13)
print('are actual logicals if 0: {}'.format((np.dot(PCX.todense(), FLAGVECS13.transpose()) % 2).sum()))

XMAT23, ZMAT23, REDQ23 = reduced_chain_complex(POSET, (2, 3))
(LX23, _), (LZ23, _) = logicals(XMAT23, ZMAT23)
print(len(LZ23))
FLAGVECS23 = np.zeros((len(LX23), NQ), dtype='uint8')
for j in range(len(LX23)):
    FLAGVECS23[j] = reduced_vec_to_flag_vec(POSET, LX23[j], REDQ23)
print('are actual logicals if 0: {}'.format((np.dot(PCX.todense(), FLAGVECS23.transpose()) % 2).sum()))

print('Checking tri-orthogonal condition for surface like X logicals: |L_j wedge L_k wedge S_l|')
REDLX = np.block([[FLAGVECS01], [FLAGVECS02], [FLAGVECS03], [FLAGVECS12], [FLAGVECS13], [FLAGVECS23]])
PCX = PCX.todense()
K, N = REDLX.shape
NSX, _ = PCX.shape
TRICOND = True
wrong = []
for ilog1, ilog2 in combinations(range(K), 2):
    for istab in range(NSX):
        TEST = np.multiply(REDLX[ilog1, :],
                           np.multiply(REDLX[ilog2, :],
                                       PCX[istab, :])).sum() % 2 == 0
        TRICOND = TRICOND and TEST
        if not TEST:
            wrong.append(((ilog1, ilog2), istab))
        if wrong:
            break
    if wrong:
        break
print('Triorthogonal: {}'.format(TRICOND))
print(REDLX[wrong[0][0][0]])
print(REDLX[wrong[0][0][1]])
print(PCX[wrong[0][1]])
