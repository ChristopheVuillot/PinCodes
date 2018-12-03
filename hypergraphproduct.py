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
                             randomclassco(checks, bits, weights))
