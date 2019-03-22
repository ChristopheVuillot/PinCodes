"""script for puncturing codes
"""
import numpy as np
import flinalg as fl
from QuantumCodeAnalysis.QuantumCodeAnalysis import low_weight_logical
from QuantumCodeAnalysis.puncturing import puncture
from QuantumCodeConstruction.utils import readsparsematrix


MX = readsparsematrix('PCMatrices/systematichp/systematic33_dim3_transpose_7_X.sms').todense()
MZ = readsparsematrix('PCMatrices/systematichp/systematic33_dim3_transpose_7_Z.sms').todense()
# MX = readsparsematrix('PCMatrices/narrowCC/narrowCC8_dim3_X.sms').todense()
# MZ = readsparsematrix('PCMatrices/narrowCC/narrowCC8_dim3_Z.sms').todense()


RX, NQ = MX.shape
RZ, _ = MZ.shape

PERM = np.random.permutation(NQ)
print('permutation: {}'.format(PERM))
K = int(NQ/60)
print('K={}'.format(K))

PUNCTMATX, PERMSTD = puncture(MX, PERM, K)
print('MX.shape = {}'.format(MX.shape))
print('PUNCTMATX.shape = {}'.format(PUNCTMATX.shape))

PUNCTSTABX = np.array(PUNCTMATX[K:, :], dtype='uint8')
PUNCTLOGX = np.array(PUNCTMATX[:K, :], dtype='uint8')
print('PUNCTSTABX.shape = {}'.format(PUNCTSTABX.shape))
print('PUNCTLOGX.shape = {}'.format(PUNCTLOGX.shape))

PUNCTSTABZ = fl.kernel(PUNCTMATX.transpose())
print('PUNCTSTABZ.shape = {}'.format(PUNCTSTABZ.shape))
KERZ = fl.kernel(PUNCTSTABX.transpose())
print('KERZ.shape = {}'.format(KERZ.shape))
PUNCTLOGZ = fl.quotient_basis(KERZ, PUNCTSTABZ)

TRIALS = 10
LOWWEIGHTLOGZ, _ = low_weight_logical(KERZ, PUNCTLOGX, TRIALS)
print('Low weight logical of weight: {}'.format(LOWWEIGHTLOGZ.sum()))
