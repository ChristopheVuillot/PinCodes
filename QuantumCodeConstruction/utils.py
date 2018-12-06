"""misc tools
"""

def writesparsematrix(matrix, filename, cardinal=2):
    """write to the file filename a binary sparse matrix
    """
    (numrows, numcolumns) = matrix.shape
    indices = zip(*matrix.nonzero())
    with open(filename, 'w') as mfile:
        mfile.write('{} {} {}\n'.format(numrows, numcolumns, cardinal))
        for row, col in indices:
            mfile.write('{} {} {}\n'.format(row+1, col+1, matrix[row, col]))
        mfile.write('0 0 0')
