"""creating pin codes from hypergraph products
"""
import numpy as np
import QuantumCodeConstruction.pincode as pinco
import QuantumCodeConstruction.hypergraphproduct as hp
from QuantumCodeConstruction.utils import writesparsematrix


if __name__ == "__main__":
    for SEED in range(1000, 1100):
        if SEED == 1097:
            # 1097 is a bad seed ! skip it !
            continue
        TRANSITIONS = hp.randomhypergraphproductlist(2, 3, 2, 2, seed=SEED, swap=True)
        POSETHP = pinco.GrPoset(TRANSITIONS, iscomplete=False)
        PCX, PCZ = pinco.pincode(POSETHP, 1, 2)
        for j, bmap in enumerate(TRANSITIONS):
            writesparsematrix(bmap, 'BoundaryMaps/randomhp_swap_{}_{}.sms'.format(SEED, j))
        WEIGHTSX = {a for che in PCX.sum(axis=1).tolist() for a in che}
        WEIGHTSZ = {a for che in PCZ.sum(axis=1).tolist() for a in che}
        print(WEIGHTSX)
        print(WEIGHTSZ)
        print(POSETHP.levelsizes)
        PCCHECKSX, PCQUBITS = PCX.shape
        PCCHECKSZ, _ = PCZ.shape
        print('qubits - xchecks - zchecks = logicals: \n \
               {} - {} - {} = {}'.format(PCQUBITS,
                                         PCCHECKSX,
                                         PCCHECKSZ,
                                         PCQUBITS - PCCHECKSX - PCCHECKSZ))

        CSSCOND = (np.dot(PCX.todense(), PCZ.transpose().todense()) % 2).sum() == 0
        print('This is a valid CSS code: {}'.format(CSSCOND))

        writesparsematrix(PCX, 'PCMatrices/randomhp_swap_{}_X.sms'.format(SEED))
        writesparsematrix(PCZ, 'PCMatrices/randomhp_swap_{}_Z.sms'.format(SEED))
