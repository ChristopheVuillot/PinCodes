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
    for _ in range(trials):
        perm = np.random.permutation(col1)
        uig1perm = uig1[:, perm]
        print(uig1perm)
        _ = fl.row_reduce_transform(uig1perm)
        print(uig1perm)
        # not col reduce no !
        redcolm = fl.col_reduce_transform(uig1perm)
        print(uig1perm)
        uig2perm = np.dot(uig2[:, perm], redcolm)
        for vec in uig1perm:
            ham_vec = sum(vec % 2)
            if ham_vec < min_d:
                if sum(np.dot(uig2perm, vec) % 2) > 0:
                    min_d = ham_vec
                    min_logical = vec
    return min_logical


def distance_upper_bound(Gx, Gz, trials):
    """
    Monte Carlo algorithm which upper bounds the distance of a CSS quantum code
    INPUT: Two linear codes given by generating matrices of shape kxn. Number of Monte Carlo trials.
    OUTPUT: Upper bound on code distance.
    """
    logicalX = low_weight_logical(Gx, Gz, trials)
    logicalZ = low_weight_logical(Gz, Gx, trials)
    dX = logicalX.hamming_weight()
    dZ = logicalZ.hamming_weight()
    return min(dX, dZ)


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
