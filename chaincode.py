""" Define posets, pinned-sets and chain codes
"""

from itertools import combinations, product
import numpy as np


class GrPoset:
    """Class for graded posets
    """

    def __init__(self, transitions, iscomplete):
        self.iscomplete = iscomplete
        self.length = len(transitions)+1
        self.transitions = [np.array(transitions[j], dtype='int') for j in range(self.length-1)]
        self.levelsizes = [np.size(transitions[j], 0) for j in range(self.length-1)] \
                          + [np.size(transitions[self.length-2], 1)]
        self.__flags__ = None
        self.__all_pinned_sets__ = [None for j in range(self.length)]
        self.__boundary_element__ = [False for j in range(self.length)]


    def neighbours_down(self, level, elem):
        """Give the neighbours below elem which is at level
        """
        assert level < self.length
        return list(np.nonzero(self.transitions[level][elem])[0])


    def neighbours_up(self, level, elem):
        """Give the neighbours above elem which is at level
        """
        assert level > 0
        return list(np.nonzero(self.transitions[level-1][:, elem])[0])


    def makecomplete(self):
        """make the poset complete by adding greatest and least elments
        """
        if not self.iscomplete:
            sumall0 = np.sum(self.transitions[0], axis=0) % 2
            if not sum(abs(sumall0)) == 0:
                self.transitions[0] = np.vstack([self.transitions[0], sumall0])
                self.levelsizes[0] += 1
                self.__boundary_element__[0] = True
            last = self.length - 2
            sumalllast = np.sum(self.transitions[last], axis=1) % 2
            if not sum(abs(sumalllast)) == 0:
                self.transitions[last] = np.transpose(np.vstack([np.transpose(self.transitions[last]), sumalllast]))
                self.levelsizes[last + 1] += 1
                self.__boundary_element__[last + 1] = True
            self.transitions = [np.ones([1, self.levelsizes[0]], dtype='int')] \
                               + self.transitions \
                               + [np.ones([self.levelsizes[last], 1], dtype='int')]
            self.levelsizes = [1] + self.levelsizes + [1]
            self.__boundary_element__ = [True] + self.__boundary_element__ + [True]
            self.length += 2
            self.iscomplete = True
            self.__all_pinned_sets__ = [None for j in range(self.length)]
            self.__flags__ = None


    def get_flags(self):
        """ get the list of the flags of the poset
        memorizes the list for latter calls
        """
        if not self.__flags__:
            flag_list = [list(range(self.levelsizes[0]))]
            for k in range(1, self.length):
                new_flag_list = []
                for flag in flag_list:
                    new_flag_list.extend([flag+[a] for a in self.neighbours_down(k-1, flag[k-1])])
                flag_list = new_flag_list
            self.__flags__ = flag_list
        return self.__flags__


    def get_all_pinned_sets(self, numberpins):
        """ get the list of the pinned sets with numberpins pins of the poset
        memorizes it for latter calls
        """
        if not self.__all_pinned_sets__[numberpins]:
            types = combinations(range(1, self.length-1), numberpins)
            pinned_sets = []
            for typ in types:
                pinss = product(*[list(range(self.levelsizes[typ[i]] - self.__boundary_element__[typ[i]])) for i in range(numberpins)])
                for pins in pinss:
                    print(pins)
                    pset = pinned_set(self.get_flags(), list(typ), list(pins))
                    if not pset == []:
                        pinned_sets.append(pset)
            self.__all_pinned_sets__[numberpins] = pinned_sets
        return self.__all_pinned_sets__[numberpins]


def projection(flag, typ):
    """project a flag according to the type typ
    """
    return [a for i, a in enumerate(flag) if i in typ]


def pinned_set(flags, typ, pins):
    """get the pinned set of type typ with pins
    """
    assert len(typ) == len(pins)
    return [f for f in flags if projection(f, typ) == pins]


def chaincode(poset, xind, zind):
    """create the parity matrices for X and Z checks
    defined by the poset and with xind and zind length
    pinned sets
    """
    poset.makecomplete()
    assert (xind + zind) < (poset.length - 2)
    flags = poset.get_flags()
    xpsets = poset.get_all_pinned_sets(xind)
    zpsets = poset.get_all_pinned_sets(zind)
    nqubit = len(flags)
    nxchecks = len(xpsets)
    nzchecks = len(zpsets)
    matx = np.zeros([nxchecks, nqubit], dtype='int')
    matz = np.zeros([nzchecks, nqubit], dtype='int')
    for checkindex, pset in enumerate(xpsets):
        for flag in pset:
            matx[checkindex][flags.index(flag)] = 1
    for checkindex, pset in enumerate(zpsets):
        for flag in pset:
            matz[checkindex][flags.index(flag)] = 1
    return (matx, matz)
