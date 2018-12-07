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
            matrix[row, col] = field(val)
        return matrix


MATX = readsparsematrix('../CCMatrices/custom_CCX.txt')
MATZ = readsparsematrix('../CCMatrices/custom_CCZ.txt')

ROWSX, QUBITSX = MATX.dimensions()
ROWSZ, QUBITSZ = MATZ.dimensions()
RANKX = MATX.rank()
RANKZ = MATZ.rank()
LOGICALS = QUBITSX - RANKX - RANKZ
MULTXZ = MATX * MATZ.transpose()
CSSCOND = sum(sum(MULTXZ)) == 0

print('mx: {} - nx:{}'.format(ROWSX, QUBITSX))
print('mx: {} - nx:{}'.format(ROWSZ, QUBITSZ))
print('rankx:{}'.format(RANKX))
print('rankz:{}'.format(RANKZ))
print('[[n = {}, k = {}]] - rate: k/n = {}'.format(QUBITSX,
                                                   LOGICALS,
                                                   float(LOGICALS)/float(QUBITSX)))
print('csscond: MX * MZ^T = 0 ? {}'.format(CSSCOND))
