"""Sage script to extract logical representatives from a CSS code
"""

import sage.all as sga

def readsparsematrix(filename):
    """Read a matrix in sparse format from file
    """
    with open(filename, 'r') as matrixfile:
        topline = matrixfile.readline()
        toplist = topline.strip().split(' ')
        numrow = int(toplist[0])
        numcol = int(toplist[1])
        fieldsize = int(toplist[2])
        field = sga.GF(fieldsize)
        matrix = sga.matrix(field, numrow, numcol, sparse=True)
        for line in matrixfile:
            linelist = line.strip().split(' ')
            row = int(linelist[0]) - 1
            col = int(linelist[1]) - 1
            val = int(linelist[2])
	    if row >= 0 and col >= 0:
	        matrix[row, col] = field(val)
        return matrix


def logicals(hmatx, hmatz):
    """Finds and returns logical operators
    from the import X and Z checks in hmatx and hmatz
    """
    _, qubitsx = hmatx.dimensions()
    ringx = hmatx.base_ring()
    assert ringx == hmatz.base_ring()
    assert qubitsx == hmatz.dimensions()[1]
    qubitvs = sga.VectorSpace(ringx, qubitsx)
    kerx = hmatx.transpose().kernel()
    imx = qubitvs.subspace(hmatx)
    kerz = hmatz.transpose().kernel()
    imz = qubitvs.subspace(hmatz)
    logicalxspace = kerz/imx
    logicalzspace = kerx/imz
    logicalx = [logicalxspace.lift(log) for log in logicalxspace.basis()]
    logicalz = [logicalzspace.lift(log) for log in logicalzspace.basis()]
    return (logicalx, logicalz)


MATX = readsparsematrix('../PCMatrices/custom_PCX.sms')
MATZ = readsparsematrix('../PCMatrices/custom_PCZ.sms')

ROWSX, QUBITSX = MATX.dimensions()
ROWSZ, QUBITSZ = MATZ.dimensions()
RANKX = MATX.rank()
RANKZ = MATZ.rank()
NUMLOGICALS = QUBITSX - RANKX - RANKZ
MULTXZ = MATX * MATZ.transpose()
CSSCOND = sum(sum(MULTXZ)) == 0

LOGX, LOGZ = logicals(MATX, MATZ)

print('mx: {} - nx:{}'.format(ROWSX, QUBITSX))
print('mx: {} - nx:{}'.format(ROWSZ, QUBITSZ))
print('rankx:{}'.format(RANKX))
print('rankz:{}'.format(RANKZ))
print('[[n = {}, k = {}]] - rate: k/n = {}'.format(QUBITSX,
                                                   NUMLOGICALS,
                                                   float(NUMLOGICALS)/float(QUBITSX)))
print('csscond: MX * MZ^T = 0 ? {}'.format(CSSCOND))
for j in range(NUMLOGICALS):
    print('logical X_{}: {}'.format(j, LOGX[j]))
print('\n')
for j in range(NUMLOGICALS):
    print('logical Z_{}: {}'.format(j, LOGZ[j]))
