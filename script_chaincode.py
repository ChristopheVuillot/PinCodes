"""testing script for chaincode.py
"""
from numpy import transpose, matmul
import chaincode as chco
import hypergraphproduct as hp

if __name__ == "__main__":
    TRANSITIONS = hp.randomhypergraphproductlist(4, 5, 2, 2, seed=None)
    POSETHP = chco.GrPoset(TRANSITIONS, iscomplete=False)
    CCX, CCZ = chco.chaincode(POSETHP, 1, 2)
    CCCHECKSX, CCQUBITS = CCX.shape
    CCCHECKSZ, _ = CCZ.shape
    print('qubits - xchecks - zchecks = logicals: \n \
            {} - {} - {} = {}'.format(CCQUBITS,
                                      CCCHECKSX,
                                      CCCHECKSZ,
                                      CCQUBITS - CCCHECKSX - CCCHECKSZ))
    print('CSS code condition: {}'.format(sum(sum(abs(matmul(CCX, transpose(CCZ)) % 2)))))
    # FLAGSHP = POSETHP.get_flags()
    # APSHP1 = POSETHP.get_all_pinned_sets(1)
    # APSHP2 = POSETHP.get_all_pinned_sets(2)
    # COUNTUP = 0
    # COUNTDOWN = 0
    # for j in range(1, POSETHP.length-1):
    #     for k in range(POSETHP.levelsizes[j]):
    #         if not POSETHP.neighbours_down(j, k):
    #             COUNTDOWN += 1
    #         if not POSETHP.neighbours_up(j, k):
    #             COUNTUP += 1
    # print(COUNTDOWN)
    # print(COUNTUP)
    # for flag in FLAGSHP:
    #     print(flag)
    # print('\n')
    # for pinset in APSHP2:
    #     print(pinset)
