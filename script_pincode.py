"""testing script for chaincode.py
"""
# import time
# import scipy.sparse as sp
import QuantumCodeConstruction.pincode as pinco
import QuantumCodeConstruction.hypergraphproduct as hp
from QuantumCodeConstruction.utils import writesparsematrix


if __name__ == "__main__":
    TRANSITIONS = hp.randomhypergraphproductlist(4, 5, 2, 2, seed=None)
    POSETHP = pinco.GrPoset(TRANSITIONS, iscomplete=False)
    PCX, PCZ = pinco.pincode(POSETHP, 1, 2)
    PCCHECKSX, PCQUBITS = PCX.shape
    PCCHECKSZ, _ = PCZ.shape
    print('qubits - xchecks - zchecks = logicals: \n \
           {} - {} - {} = {}'.format(PCQUBITS,
                                     PCCHECKSX,
                                     PCCHECKSZ,
                                     PCQUBITS - PCCHECKSX - PCCHECKSZ))
    # TIME = time.localtime()
    writesparsematrix(PCX, 'PCMatrices/custom_PCX.sms')
    writesparsematrix(PCZ, 'PCMatrices/custom_PCZ.sms')
