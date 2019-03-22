"""functions to puncture codes
"""

import numpy as np
import flinalg as fl


def puncture(hmat, perm, k):
    """puncture the code given by hmat
    after permuting the columns according to perm
    and reducing it to standard form.
    return the punctured matrix as well as the additional
    permutation done during standard form process.
    """
    hmatcopy = np.array(hmat, dtype='uint8')[:, perm]
    _, permstd = fl.standard_form(hmatcopy)
    return hmatcopy[:, k:], permstd
