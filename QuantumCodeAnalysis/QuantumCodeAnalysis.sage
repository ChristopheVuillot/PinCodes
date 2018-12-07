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

MX = readsparsematrix('../CCMatrices/custom_CCX.txt')
MZ = readsparsematrix('../CCMatrices/custom_CCZ.txt')

mx, nx = MX.dimensions()
mz, nz = MZ.dimensions()
rankx = MX.rank()
rankz = MZ.rank()
k = nx - rankx - rankz
MultXZ = MX * MZ.transpose()
csscond = sum(sum(MultXZ)) == 0

print('mx: {} - nx:{}'.format(mx, nx))
print('mx: {} - nx:{}'.format(mz, nz))
print('rankx:{}'.format(rankx))
print('rankz:{}'.format(rankz))
print('[[n = {}, k = {}]] - rate: k/n = {}'.format(nx, k, float(k)/float(nx)))
print('csscond: MX * MZ^T = 0 ? {}'.format(csscond))
