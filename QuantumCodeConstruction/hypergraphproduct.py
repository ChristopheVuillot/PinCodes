"""Constructing random hypergraph product codes
"""

import numpy as np
import numpy.random as rd


def systematicclassco(checks, bits):
    """returns the list of all, 2^(checks*bits), checks x bits binary matrices.
    removing emty rows and empty columns.
    use with care.
    """
    nbits = checks*bits
    nmatrices = 2**nbits
    hmatlist = [np.zeros((checks, bits), dtype='uint8') for _ in range(1, nmatrices)]
    for j in range(1, nmatrices):
        bitstring = np.binary_repr(j, width=nbits)
        for row in range(checks):
            for col in range(bits):
                hmatlist[j-1][row, col] = int(bitstring[row*bits + col])
        emptyrows = []
        emptycolumns = []
        for k in range(checks):
            if not hmatlist[j-1][k].any():
                emptyrows.append(k)
        for k in range(bits):
            if not hmatlist[j-1][:, k].any():
                emptycolumns.append(k)
        hmatlist[j-1] = np.delete(hmatlist[j-1], emptyrows, 0)
        hmatlist[j-1] = np.delete(hmatlist[j-1], emptycolumns, 1)
    return hmatlist


def randomclassco(checks, bits, weights):
    """construct randomly a c by n adjacency matrix
    with an average of w ones per line
    """
    hmat = np.zeros((checks, bits), dtype='int')
    for i in range(0, checks):
        for j in range(0, bits):
            if rd.random() < weights/bits:
                hmat[i, j] = 1
    emptyrows = []
    emptycolumns = []
    for j in range(checks):
        if not hmat[j].any():
            emptyrows.append(j)
    for j in range(bits):
        if not hmat[:, j].any():
            emptycolumns.append(j)
    hmat = np.delete(hmat, emptyrows, 0)
    hmat = np.delete(hmat, emptycolumns, 1)
    return hmat


def hypergraphproduct(hmat1, hmat2):
    """Construct the hypergraph product from
    two given classical codes
    """
    (checks1, bits1) = hmat1.shape
    (checks2, bits2) = hmat2.shape
    hmatx = np.bmat([np.kron(np.identity(checks2), hmat1),
                     np.kron(hmat2, np.identity(checks1))])
    hmatz = np.bmat([np.kron(np.transpose(hmat2), np.identity(bits1)),
                     np.kron(-np.identity(bits2), np.transpose(hmat1))])
    return (hmatx, hmatz)


def randomhypergraphproduct(checks, bits, weights, seed=None):
    """create a random quantum code
    from two random classical codes
    """
    rd.seed(seed)
    return hypergraphproduct(randomclassco(checks, bits, weights),
                             np.transpose(randomclassco(checks, bits, weights)))


def hypergraphproductlist(hmatlist, hmatnew):
    """Construct the product from a chain complex and
    an additional classical code
    """
    length = len(hmatlist)
    checkbitlist = [hmat.shape for hmat in hmatlist]
    (checks, bits) = hmatnew.shape
    transitionlist = [np.bmat([np.kron(np.identity(checks), hmatlist[0]),
                               np.kron(hmatnew, np.identity(checkbitlist[0][0]))])]
    for j in range(1, length):
        transitionlist.append(np.bmat([[np.kron(np.identity(checks), hmatlist[j]),
                                        np.kron(hmatnew, np.identity(checkbitlist[j][0]))],
                                       [np.zeros((bits*checkbitlist[j-1][0], checks*checkbitlist[j][1]), dtype='int'),
                                        np.kron(np.identity(bits), hmatlist[j-1])]]))

    transitionlist.append(np.bmat([[np.kron(hmatnew, np.identity(checkbitlist[length-1][1]))],
                                   [np.kron(np.identity(bits), hmatlist[length-1])]]))
    return transitionlist


def reapeatedhypergraphproduct(hmat, nrepeat, transpose=True):
    """produce the repeated hypergraph product of the matrice hmat
    alternatively transposing it if transpose.
    """
    transitionlist = [hmat]
    for j in range(nrepeat):
        if transpose and ((j % 2) == 0):
            transitionlist = hypergraphproductlist(transitionlist, hmat.transpose())
        else:
            transitionlist = hypergraphproductlist(transitionlist, hmat)
    return transitionlist


def hypergraphproductlistlist(hmatlist):
    """produce the hypergraph product of all elements in hmatlist
    """
    transitionlist = [hmatlist[0]]
    for j in range(1, len(hmatlist)):
        transitionlist = hypergraphproductlist(transitionlist, hmatlist[j])
    return transitionlist


def randomhypergraphproductlist(checks, bits, weights, length, seed=None, swap=False):
    """Iterates the product construction with random classical codes
    length times
    """
    rd.seed(seed)
    transitions = [randomclassco(checks, bits, weights)]
    for _ in range(length):
        if swap:
            (checks, bits) = (bits, checks)
        transitions = hypergraphproductlist(transitions, randomclassco(checks, bits, weights))
    return transitions
