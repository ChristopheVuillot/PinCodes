"""Sage functions to extract logical representatives from a CSS code
"""

import sage.all as sga
from itertools import combinations
from scipy.sparse import load_npz


def readsparsematrix(filename, fieldsize=2):
    """Read a matrix in sparse format from file
    """
    _, ext = filename.split('.')
    if ext == 'sms':
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
    elif ext == 'npz':
        field = sga.GF(fieldsize)
        matrix = sga.matrix(field, load_npz(filename).todense())
    else:
        print('Wrong format !! .{} is not recognized'.format(ext))
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
    logicalxspace = kerz / imx
    logicalzspace = kerx / imz
    logicalx = [logicalxspace.lift(log) for log in logicalxspace.basis()]
    logicalz = [logicalzspace.lift(log) for log in logicalzspace.basis()]
    return (logicalx, logicalz)


def logical_circuit(logicalx, k_orth):
    """find the kth level equivalent logical circuit
    for a transversal application of R(k_orth) on
    the code given all the logical X operators
    """
    k = len(logicalx)
    logicalx_index = zip(range(k), logicalx)
    gates = [[] for _ in range(k_orth)]
    for j, jgates in enumerate(gates):
        for tuple in combinations(logicalx_index, j + 1):
            indices, logicals = zip(*tuple)
            coef = ((-1)**j * 2**j * sum([int(prod(t)) for t in zip(*logicals)])) \
                % (2**k_orth)
            if coef != 0:
                jgates.append((indices, coef))
    return gates
