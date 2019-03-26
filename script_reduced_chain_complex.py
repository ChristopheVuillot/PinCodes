"""script to test reduced chain complexes of pin codes
"""
import numpy as np
from QuantumCodeConstruction.pincode import GrPoset, reduced_chain_complex, pincode
from QuantumCodeConstruction.utils import readsparsematrix
from QuantumCodeAnalysis.QuantumCodeAnalysis import low_weight_logical, logicals


# TRANSITIONS = [readsparsematrix('BoundaryMaps/systematichp/systematic33_dim3_transpose_7_{}.sms'.format(j)).todense() for j in range(3)]
TRANSITIONS = [readsparsematrix('BoundaryMaps/535_6840_{}.npz'.format(j)).todense().transpose() for j in range(1, 4)]
POSET = GrPoset(TRANSITIONS, iscomplete=False)
PCX, PCZ = pincode(POSET, 1, 2)
(LX, _), (LZ, _) = logicals(PCX.todense(), PCZ.todense())
print(len(LX))
print(PCX.todense())
print(PCZ.todense())
XMAT0, ZMAT0, REDQ0 = reduced_chain_complex(POSET, (0,))
(LX0, _), (LZ0, _) = logicals(XMAT0, ZMAT0)
XMAT1, ZMAT1, REDQ1 = reduced_chain_complex(POSET, (1,))
(LX1, _), (LZ1, _) = logicals(XMAT1, ZMAT1)
XMAT2, ZMAT2, REDQ2 = reduced_chain_complex(POSET, (2,))
(LX2, _), (LZ2, _) = logicals(XMAT2, ZMAT2)
XMAT3, ZMAT3, REDQ3 = reduced_chain_complex(POSET, (3,))
(LX3, _), (LZ3, _) = logicals(XMAT3, ZMAT3)
print(len(LX0))
print(len(LX1))
print(len(LX2))
print(len(LX3))
XMAT01, ZMAT01, REDQ01 = reduced_chain_complex(POSET, (0, 1))
(LX01, _), (LZ01, _) = logicals(XMAT01, ZMAT01)
XMAT02, ZMAT02, REDQ02 = reduced_chain_complex(POSET, (0, 2))
(LX02, _), (LZ02, _) = logicals(XMAT02, ZMAT02)
XMAT03, ZMAT03, REDQ03 = reduced_chain_complex(POSET, (0, 3))
(LX03, _), (LZ03, _) = logicals(XMAT03, ZMAT03)
XMAT12, ZMAT12, REDQ12 = reduced_chain_complex(POSET, (1, 2))
(LX12, _), (LZ12, _) = logicals(XMAT12, ZMAT12)
XMAT13, ZMAT13, REDQ13 = reduced_chain_complex(POSET, (1, 3))
(LX13, _), (LZ13, _) = logicals(XMAT13, ZMAT13)
XMAT23, ZMAT23, REDQ23 = reduced_chain_complex(POSET, (2, 3))
(LX23, _), (LZ23, _) = logicals(XMAT23, ZMAT23)
print(len(LZ01))
print(len(LZ02))
print(len(LZ03))
print(len(LZ12))
print(len(LZ13))
print(len(LZ23))
