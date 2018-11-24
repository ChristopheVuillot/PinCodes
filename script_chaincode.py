"""testing script for chaincode.py
"""
from numpy import transpose, matmul
from numpy.linalg import matrix_rank
import chaincode as chco
import hypergraphproduct as hp

# TRAN = [[[1, 1, 1]],
#         [[1, 1, 0], [1, 0, 1], [0, 1, 1]],
#         [[1], [1], [1]]]
# POSET = chco.GrPoset(TRAN, True)
#
# Graphical representation
#           0
#          /|\
#         0 1 2
#         |X X|
#         0 1 2
#          \|/
#           0
#
# FLAGS = POSET.get_flags()
# APS = POSET.get_all_pinned_sets(2)
# print(POSET.levelsizes)
# print(FLAGS)
# print([chco.projection(f, [0, 2]) for f in FLAGS])
# print(chco.pinned_set(FLAGS, [1], [2]))
# print([[FLAGS.index(f) for f in ps] for ps in APS])
# print(POSET.length)
# MX, MZ = chco.chaincode(POSET, 1, 0)
# print(matmul(MX, transpose(MZ)) % 2)

HPMX, HPMZ = hp.randomhypergraphproduct(10, 12, 4)
CHECKSX, QUBITS = HPMX.shape
CHECKSZ, _ = HPMZ.shape
POSETHP = chco.GrPoset([HPMX, transpose(HPMZ)], iscomplete=False)

CCX, CCZ = chco.chaincode(POSETHP, 1, 1)
FLAGSHP = POSETHP.get_flags()
APSHP = POSETHP.get_all_pinned_sets(1)
print(CCX.shape)
print(matrix_rank(CCX))
CCCHECKS, CCQUBITS = CCX.shape
print(sum(sum(abs(matmul(CCX, transpose(CCZ)) % 2))))
