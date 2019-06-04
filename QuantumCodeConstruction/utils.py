"""misc tools
"""

from scipy.sparse import load_npz, save_npz, lil_matrix


def writesparsematrix(matrix, filename, cardinal=2):
    """write to the file filename a binary sparse matrix
    """
    _, ext = filename.split('.')
    if ext == 'sms':
        (numrows, numcolumns) = matrix.shape
        indices = zip(*matrix.nonzero())
        with open(filename, 'w') as mfile:
            mfile.write('{} {} {}\n'.format(numrows, numcolumns, cardinal))
            for row, col in indices:
                mfile.write('{} {} {}\n'.format(row + 1, col + 1, matrix[row, col]))
            mfile.write('0 0 0')
    elif ext == 'npz':
        save_npz(matrix.tocsc(), filename)
    else:
        print('Wrong format !! .{} is not recognized'.format(ext))


def readsparsematrix(filename):
    """Read a matrix in sparse format from file
    """
    _, ext = filename.split('.')
    if ext == 'sms':
        with open(filename, 'r') as matrixfile:
            topline = matrixfile.readline()
            toplist = topline.strip().split(' ')
            numrow = int(toplist[0])
            numcol = int(toplist[1])
            matrix = lil_matrix((numrow, numcol), dtype='uint8')
            for line in matrixfile:
                linelist = line.strip().split(' ')
                row = int(linelist[0]) - 1
                col = int(linelist[1]) - 1
                val = int(float(linelist[2]))
                if row >= 0 and col >= 0:
                    matrix[row, col] += val
    elif ext == 'npz':
        matrix = load_npz(filename)
    else:
        print('Wrong format !! .{} is not recognized'.format(ext))
    return matrix
