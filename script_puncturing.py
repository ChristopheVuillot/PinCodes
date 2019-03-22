"""script for puncturing codes
"""
import numpy as np
import flinalg as fl
from QuantumCodeAnalysis.QuantumCodeAnalysis import low_weight_logical
from QuantumCodeAnalysis.puncturing import puncture
from QuantumCodeConstruction.utils import readsparsematrix


# MX = readsparsematrix('PCMatrices/systematichp/systematic33_dim3_transpose_429_X.sms').todense()
# MZ = readsparsematrix('PCMatrices/systematichp/systematic33_dim3_transpose_429_Z.sms').todense()
MX = readsparsematrix('PCMatrices/narrowCC/narrowCC2_dim6_X.sms').todense()
MZ = readsparsematrix('PCMatrices/narrowCC/narrowCC2_dim6_Z.sms').todense()
# MX = readsparsematrix('PCMatrices/narrowCC/narrowCC_alt42_dim6_X.sms').todense()
# MZ = readsparsematrix('PCMatrices/narrowCC/narrowCC_alt42_dim6_Z.sms').todense()
# MX = readsparsematrix('PCMatrices/535_3420_XCOS.sms').todense()
# MZ = readsparsematrix('PCMatrices/535_3420_ZCOS.sms').todense()


RX, NQ = MX.shape
RZ, _ = MZ.shape

PERM = np.random.permutation(NQ)
print('permutation: {}'.format(PERM))
K = int(NQ/10)
print('K={}'.format(K))

PUNCTMATX, PERMSTD = puncture(MX, PERM, K)
print(PERMSTD)
print('MX.shape = {}'.format(MX.shape))
print('PUNCTMATX.shape = {}'.format(PUNCTMATX.shape))

PUNCTSTABX = np.array(PUNCTMATX[K:, :], dtype='uint8')
PUNCTLOGX = np.array(PUNCTMATX[:K, :], dtype='uint8')
print('PUNCTSTABX.shape = {}'.format(PUNCTSTABX.shape))
print('PUNCTLOGX.shape = {}'.format(PUNCTLOGX.shape))

PUNCTSTABZ = fl.kernel(PUNCTMATX.transpose())
print('PUNCTSTABZ.shape = {}'.format(PUNCTSTABZ.shape))
KERZ = fl.kernel(PUNCTSTABX.transpose())
KERX = fl.kernel(PUNCTSTABZ.transpose())

print('KERZ.shape = {}'.format(KERZ.shape))
print('KERX.shape = {}'.format(KERX.shape))
PUNCTLOGZ = fl.quotient_basis(KERZ, PUNCTSTABZ)
REALLOGX = fl.quotient_basis(KERX, PUNCTSTABX)
print('REALLOGX len = {}'.format(len(REALLOGX)))
print('PUNCTLOGZ len = {}'.format(len(PUNCTLOGZ)))

TRIALS = 10
LOWWEIGHTLOGZ, _ = low_weight_logical(KERZ, PUNCTLOGX, TRIALS)
print('Low weight logical of weight: {}'.format(LOWWEIGHTLOGZ.sum()))
