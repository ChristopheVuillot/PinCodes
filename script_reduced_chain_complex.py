"""script to test reduced chain complexes of pin codes
"""
import numpy as np
from QuantumCodeConstruction.pincode import GrPoset, reduced_chain_complex, pincode
from QuantumCodeConstruction.utils import readsparsematrix
from QuantumCodeAnalysis.QuantumCodeAnalysis import low_weight_logical, logicals


TRANSITIONS = [readsparsematrix('BoundaryMaps/systematichp/systematic33_dim3_transpose_1_{}.sms'.format(j)).todense() for j in range(3)]
POSET = GrPoset(TRANSITIONS, iscomplete=False)
PCX, PCZ = pincode(POSET, 1, 2)
XMAT2, ZMAT2, REDQ2 = reduced_chain_complex(POSET, (2,))
(LX2, _), (LZ2, _) = logicals(XMAT2, ZMAT2)
print(LX2)
print(LZ2)
print(XMAT2)
print(ZMAT2)
print(REDQ2)
XMAT12, ZMAT12, REDQ12 = reduced_chain_complex(POSET, (1, 2))
(LX12, _), (LZ12, _) = logicals(XMAT12, ZMAT12)
print(LX12)
print(LZ12)
print(XMAT12)
print(ZMAT12)
print(REDQ12)
