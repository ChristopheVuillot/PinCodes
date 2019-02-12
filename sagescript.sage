load('QuantumCodeAnalysis/QuantumCodeAnalysis.sage')

MATX = readsparsematrix('PCMatrices/535_6840_X.sms')
MATZ = readsparsematrix('PCMatrices/535_6840_Z.sms')

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
