import numpy as np
from sage.all import GF, matrix, vector, Permuations, copy, LinearCode


def low_weight_logical(G1, G2, trials):
    """
    Monte Carlo algorithm which tries to find the lowest weight code word of G1 which has odd overlap with
    code words in G2.
    INPUT: Two linear codes given by generating matrices of shape kxn.
    OUTPUT: A code word of G1 with low weight.
    """
    k, n = G1.dimensions()
    min_d = n
    min_logical = None
    for i in range(trials):
        perm_rand = Permutations(n).random_element()
        G1perm = copy(G1)
        G1perm.permute_columns(perm_rand)
        C1 = LinearCode(G1perm)
        C1, perm_stdform = C1.standard_form(return_permutation=True)
        G1perm = copy(C1.systematic_generator_matrix())
        G1perm.permute_columns((perm_rand*perm_stdform).inverse())
        for v in G1perm:
            ham_v = v.hamming_weight()
            if ham_v < min_d:
                if (G2 * v).hamming_weight() > 0:
                    min_d = ham_v
                    min_logical = v
    return min_logical


def distance_upper_bound(Gx, Gz, trials):
    logicalX = low_weight_logical(Gx, Gz, trials)
    logicalZ = low_weight_logical(Gz, Gx, trials)
    dX = logicalX.hamming_weight()
    dZ = logicalZ.hamming_weight()
    return min(dX, dZ)
