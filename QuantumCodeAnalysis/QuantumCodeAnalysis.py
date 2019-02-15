"""Functions to analyse quantum codes
"""
from itertools import combinations
import numpy as np
# from sage.all import GF, matrix, vector, Permuations, copy, LinearCode
import flinalg as fl

def low_weight_logical(gen1, gen2, trials):
    """
    Monte Carlo algorithm which tries to find the lowest weight code word of G1
    which has odd overlap with code words in G2.
    INPUT: Two linear codes given by generating matrices of shape kxn. Number of Monte Carlo trials.
    OUTPUT: A code word of G1 with low weight.
    """
    uig1 = np.array(gen1, dtype='uint8')
    uig2 = np.array(gen2, dtype='uint8')
    _, col1 = uig1.shape
    min_d = col1
    min_logical = None
    total_permutation = np.arange(col1)
    for _ in range(trials):
        perm = np.random.permutation(col1)
        uig1 = uig1[:, perm]
        _, permstdf = fl.standard_form(uig1)
        uig2 = uig2[:, perm][:, permstdf]
        total_permutation = total_permutation[perm][permstdf]
        for vec in uig1:
            ham_vec = sum(vec % 2)
            if ham_vec < min_d:
                if sum(np.dot(uig2, vec) % 2) > 0:
                    min_d = ham_vec
                    min_logical = vec
    return min_logical, total_permutation


def distance_upper_bound(logicalx, checkx, logicalz, checkz, trials):
    """
    Monte Carlo algorithm which upper bounds the distance of a CSS quantum code
    INPUT: Two linear codes given by generating matrices of shape kxn. Number of Monte Carlo trials.
    OUTPUT: Upper bound on code distance.
    """
    lowx, _ = low_weight_logical(np.block([[logicalx], [checkx]]), logicalz, trials)
    lowz, _ = low_weight_logical(np.block([[logicalz], [checkz]]), logicalx, trials)
    disx = lowx.hamming_weight()
    disz = lowz.hamming_weight()
    return min(disx, disz)


def logicals(hmatx, hmatz):
    """Finds and returns logical operators
    from the import X and Z checks in hmatx and hmatz
    """
    uintmatx = np.array(hmatx, dtype='uint8')
    uintmatz = np.array(hmatz, dtype='uint8')
    kerx = fl.kernel(np.transpose(uintmatx))
    kerz = fl.kernel(np.transpose(uintmatz))
    logicalxspace = fl.quotient_basis(kerz, uintmatx)
    logicalzspace = fl.quotient_basis(kerx, uintmatz)
    return (logicalxspace, logicalzspace)


def logical_circuit(logicalx, k_orth):
    """find the kth level equivalent logical circuit
    for a transversal application of R(k_orth) on
    the code given all the logical X operators
    """
    k = len(logicalx)
    logicalx_index = list(zip(range(k), logicalx))
    gates = [[] for _ in range(k_orth)]
    for j, jgates in enumerate(gates):
        print(j)
        for tupl in combinations(logicalx_index, j + 1):
            indices, logs = zip(*tupl)
            print(indices)
            wedge = [int(np.prod(t)) for t in zip(*logs)]
            print(wedge)
            coef = ((-1)**j * 2**j * sum(wedge)) \
                % (2**k_orth)
            print(coef)
            if coef != 0:
                jgates.append((indices, coef))
    return gates
