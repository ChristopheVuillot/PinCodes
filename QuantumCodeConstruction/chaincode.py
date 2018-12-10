""" Define posets, pinned-sets and chain codes
"""

from itertools import combinations, product
import scipy.sparse as sp
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
        self.__all_pinned_sets__ = [None for _ in range(self.length)]
        self.__boundary_element__ = [False for _ in range(self.length)]
        self.__all_neighbours_down__ = [[None for _ in range(self.levelsizes[j])] for j in range(self.length)]
        self.__all_neighbours_up__ = [[None for _ in range(self.levelsizes[j])] for j in range(self.length)]


    def neighbours_down(self, level, elem):
        """Give the neighbours below elem which is at level
        """
        assert level < self.length
        assert elem < self.levelsizes[level]
        if not self.__all_neighbours_down__[level][elem]:
            self.__all_neighbours_down__[level][elem] = list(np.nonzero(self.transitions[level][elem])[0])
        return self.__all_neighbours_down__[level][elem]


    def neighbours_up(self, level, elem):
        """Give the neighbours above elem which is at level
        """
        assert level > 0
        assert elem < self.levelsizes[level]
        if not self.__all_neighbours_up__[level][elem]:
            self.__all_neighbours_up__[level][elem] = list(np.nonzero(self.transitions[level-1][:, elem])[0])
        return self.__all_neighbours_up__[level][elem]


    def makecomplete(self):
        """make the poset complete by adding greatest and least elments
        """
        # TODO:
        # check if there are hanging elements in the middle and patch them
        # to link them to level 0 or D.
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
            self.__all_neighbours_down__ = [[None for _ in range(self.levelsizes[j])]
                                            for j in range(self.length)]
            self.__all_neighbours_up__ = [[None for _ in range(self.levelsizes[j])]
                                          for j in range(self.length)]


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
                    pset = self.pinned_set(list(typ), list(pins))
                    # pset = pinned_set_simple(self.get_flags(), list(typ), list(pins))
                    if not pset == []:
                        pinned_sets.append(pset)
            self.__all_pinned_sets__[numberpins] = pinned_sets
        return self.__all_pinned_sets__[numberpins]


    def pinned_set(self, typ, pins):
        """get the pinned set of type typ with pins
        """
        numpins = len(pins)
        # assert len(typ) == numpins
        # assert numpins < self.length
        pset = [[pins[0]]]
        for level in range(typ[0], 0, -1):
            new_pset = []
            for flag in pset:
                new_pset.extend([[a]+flag for a in self.neighbours_up(level, flag[0])])
            pset = new_pset
        for interval in range(numpins-1):
            for level in range(typ[interval], typ[interval+1]-1):
                new_pset = []
                for flag in pset:
                    new_pset.extend([flag+[a] for a in self.neighbours_down(level, flag[level])])
                pset = new_pset
            pset = [flag+[pins[interval+1]] for flag in pset
                    if pins[interval+1] in self.neighbours_down(typ[interval+1]-1,
                                                                flag[typ[interval+1]-1])]
            if not pset:
                return []
        for level in range(typ[numpins-1], self.length-1):
            new_pset = []
            for flag in pset:
                new_pset.extend([flag+[a] for a in self.neighbours_down(level, flag[level])])
            pset = new_pset
        return pset


def projection(flag, typ):
    """project a flag according to the type typ
    """
    return [a for i, a in enumerate(flag) if i in typ]


def pinned_set_simple(flags, typ, pins):
    """get the pinned set of type typ with pins
    """
    # assert len(typ) == len(pins)
    return [f for f in flags if projection(f, typ) == pins]


def chaincode(poset, xind, zind):
    """create the parity matrices for X and Z checks
    defined by the poset and with xind and zind length
    pinned sets
    """
    poset.makecomplete()
    assert (xind + zind) < (poset.length - 2)
    flags = poset.get_flags()
    # print(flags)
    xpsets = poset.get_all_pinned_sets(xind)
    zpsets = poset.get_all_pinned_sets(zind)
    nqubit = len(flags)
    nxchecks = len(xpsets)
    nzchecks = len(zpsets)
    matx = sp.dok_matrix((nxchecks, nqubit), dtype='int')
    matz = sp.dok_matrix((nzchecks, nqubit), dtype='int')
    for checkindex, pset in enumerate(xpsets):
        for flag in pset:
            matx[checkindex, flags.index(flag)] = 1
    for checkindex, pset in enumerate(zpsets):
        for flag in pset:
            matz[checkindex, flags.index(flag)] = 1
    return (matx, matz)
