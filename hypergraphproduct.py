"""Constructing random hypergraph product codes
"""

import numpy as np
import numpy.random as rd


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
                     -np.kron(np.identity(bits2), np.transpose(hmat1))])
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


def randomhypergraphproductlist(checks, bits, weights, length, seed=None):
    """Iterates the product construction with random classical codes
    length times
    """
    rd.seed(seed)
    transitions = [randomclassco(checks, bits, weights)]
    for _ in range(length):
        transitions = hypergraphproductlist(transitions, randomclassco(checks, bits, weights))
    return transitions
