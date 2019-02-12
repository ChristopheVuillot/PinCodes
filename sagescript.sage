load('QuantumCodeAnalysis/QuantumCodeAnalysis.sage')

MATX = readsparsematrix('PCMatrices/535_3420_X.sms')
MATZ = readsparsematrix('PCMatrices/535_3420_Z.sms')
MATXCOS = readsparsematrix('PCMatrices/535_3420_XCOS.sms')
MATZCOS = readsparsematrix('PCMatrices/535_3420_ZCOS.sms')

ROWSX, QUBITSX = MATX.dimensions()
ROWSZ, QUBITSZ = MATZ.dimensions()
RANKX = MATX.rank()
RANKZ = MATZ.rank()
NUMLOGICALS = QUBITSX - RANKX - RANKZ
MULTXZ = MATX * MATZ.transpose()
CSSCOND = sum(sum(matrix(ZZ, MULTXZ))) == 0

print('mx: {} - nx:{}'.format(ROWSX, QUBITSX))
print('mz: {} - nz:{}'.format(ROWSZ, QUBITSZ))
print('rankx:{}'.format(RANKX))
print('rankz:{}'.format(RANKZ))
print('[[n = {}, k = {}]]'.format(QUBITSX, NUMLOGICALS))
print('rate: k/n = {}'.format(float(NUMLOGICALS) / float(QUBITSX)))
print('csscond: MX * MZ^T = 0 ? {}'.format(CSSCOND))

ROWSXCOS, QUBITSXCOS = MATXCOS.dimensions()
ROWSZCOS, QUBITSZCOS = MATZCOS.dimensions()
RANKXCOS = MATXCOS.rank()
RANKZCOS = MATZCOS.rank()
NUMLOGICALSCOS = QUBITSXCOS - RANKXCOS - RANKZCOS
MULTXZCOS = MATXCOS * MATZCOS.transpose()
CSSCONDCOS = sum(sum(matrix(ZZ, MULTXZCOS))) == 0

print('mx: {} - nx:{}'.format(ROWSXCOS, QUBITSXCOS))
print('mz: {} - nz:{}'.format(ROWSZCOS, QUBITSZCOS))
print('rankx:{}'.format(RANKXCOS))
print('rankz:{}'.format(RANKZCOS))
print('[[n = {}, k = {}]]'.format(QUBITSXCOS, NUMLOGICALSCOS))
print('rate: k/n = {}'.format(float(NUMLOGICALSCOS) / float(QUBITSXCOS)))
print('csscond: MX * MZ^T = 0 ? {}'.format(CSSCONDCOS))


LOGX, LOGZ = logicals(MATX, MATZ)

LOGCIRC = logical_circuit(LOGX, 3)
# print(LOGCIRC)
for gate in LOGCIRC[0]:
    if 1 in gate[0]:
        print(gate)
for gate in LOGCIRC[1]:
    if 1 in gate[0]:
        print(gate)
for gate in LOGCIRC[2]:
    if 1 in gate[0]:
        print(gate)
# for j in range(NUMLOGICALS):
#     print('logical X_{}: {}'.format(j, LOGX[j]))
# print('\n')
# for j in range(NUMLOGICALS):
#     print('logical Z_{}: {}'.format(j, LOGZ[j]))
