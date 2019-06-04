"""functions to unpuncture codes
"""

import numpy as np
import flinalg as fl


def unpuncture(hmat, LOGX):
    """unpuncture the code given by hmat
    by adding X-logicals
    """
    k = LOGX.shape[0]
    r = hmat.shape[0]
    return np.block([
                     [np.eye(k,dtype='uint8'), LOGX],
                     [np.zeros((r,k),dtype='uint8'), hmat]
                    ])
