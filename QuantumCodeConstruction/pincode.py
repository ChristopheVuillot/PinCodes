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
        self.length = len(transitions) + 1
        self.transitions = [np.array(transitions[j], dtype='int')
                            for j in range(self.length - 1)]
        self.levelsizes = [np.size(transitions[j], 0)
                           for j in range(self.length - 1)] + [np.size(transitions[self.length - 2], 1)]
        self.__flags__ = None
        self.__all_pinned_sets__ = [None for _ in range(self.length)]
        self.__all_pinned_sets_with_bound__ = [None for _ in range(self.length)]
        self.__boundary_element__ = [False for _ in range(self.length)]
        self.__all_neighbours_down__ = [[None for _ in range(self.levelsizes[j])]
                                        for j in range(self.length)]
        self.__all_neighbours_up__ = [[None for _ in range(self.levelsizes[j])]
                                      for j in range(self.length)]

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
            self.__all_neighbours_up__[level][elem] = list(np.nonzero(self.transitions[level - 1][:, elem])[0])
        return self.__all_neighbours_up__[level][elem]

    def makecomplete(self):
        """make the poset complete by adding greatest and least elments
        """
        if not self.iscomplete:
            sumall0 = np.sum(self.transitions[0], axis=0) % 2
            if not sum(abs(sumall0)) == 0:
                print("Adding boundary element !")
                self.transitions[0] = np.vstack([self.transitions[0], sumall0])
                self.levelsizes[0] += 1
                self.__boundary_element__[0] = True
            last = self.length - 2
            sumalllast = np.sum(self.transitions[last], axis=1) % 2
            if not sum(abs(sumalllast)) == 0:
                print("Adding boundary element !")
                self.transitions[last] = np.transpose(np.vstack([np.transpose(self.transitions[last]), sumalllast]))
                self.levelsizes[last + 1] += 1
                self.__boundary_element__[last + 1] = True
            self.iscomplete = True
            self.__all_pinned_sets__ = [None for j in range(self.length)]
            self.__all_pinned_sets_with_bound__ = [None for j in range(self.length)]
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
            flag_list = [[a] for a in range(self.levelsizes[0])]
            for k in range(1, self.length):
                new_flag_list = []
                for flag in flag_list:
                    new_flag_list.extend([flag + [a] for a in self.neighbours_down(k - 1, flag[k - 1])])
                flag_list = new_flag_list
            self.__flags__ = flag_list
        return self.__flags__

    def get_all_pinned_sets(self, numberpins):
        """ get the list of the pinned sets with numberpins pins of the poset
        memorizes it for latter calls
        """
        if not self.__all_pinned_sets__[numberpins]:
            types = combinations(range(0, self.length), numberpins)
            pinned_sets = {}
            for typ in types:
                pinss = product(*[list(range(self.levelsizes[typ[i]] - self.__boundary_element__[typ[i]])) for i in range(numberpins)])
                for pins in pinss:
                    pset = self.pinned_set(list(typ), list(pins))
                    if not pset == []:
                        pinned_sets.setdefault(typ, []).append(pset)
            self.__all_pinned_sets__[numberpins] = pinned_sets
        return self.__all_pinned_sets__[numberpins]

    def get_all_pinned_sets_with_bound(self, numberpins):
        """ get the list of the pinned sets with numberpins pins of the poset
        memorizes it for latter calls
        """
        if not self.__all_pinned_sets_with_bound__[numberpins]:
            types = combinations(range(0, self.length), numberpins)
            pinned_sets = {}
            for typ in types:
                pinss = product(*[list(range(self.levelsizes[typ[i]])) for i in range(numberpins)])
                for pins in pinss:
                    pset = self.pinned_set(list(typ), list(pins))
                    if not pset == []:
                        pinned_sets.setdefault(typ, []).append(pset)
            self.__all_pinned_sets_with_bound__[numberpins] = pinned_sets
        return self.__all_pinned_sets_with_bound__[numberpins]

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
                new_pset.extend([[a] + flag for a in self.neighbours_up(level, flag[0])])
            pset = new_pset
        for interval in range(numpins - 1):
            for level in range(typ[interval], typ[interval + 1] - 1):
                new_pset = []
                for flag in pset:
                    new_pset.extend([flag + [a] for a in self.neighbours_down(level, flag[level])])
                pset = new_pset
            pset = [flag + [pins[interval + 1]] for flag in pset
                    if pins[interval + 1] in self.neighbours_down(typ[interval + 1] - 1,
                                                                  flag[typ[interval + 1] - 1])]
            if not pset:
                return []
        for level in range(typ[numpins - 1], self.length - 1):
            new_pset = []
            for flag in pset:
                new_pset.extend([flag + [a] for a in self.neighbours_down(level, flag[level])])
            pset = new_pset
        return pset

    def check_boundary_map(self):
        """checks that the boundary map is correct,
        i.e it is zero when squared.
        """
        isboundmap = True
        for j in range(self.length - 2):
            if (np.dot(self.transitions[j], self.transitions[j + 1]).sum() % 2) != 0:
                isboundmap = False
                break
        return isboundmap


def projection(flag, typ):
    """project a flag according to the type typ
    """
    return [a for i, a in enumerate(flag) if i in typ]


def pinned_set_simple(flags, typ, pins):
    """get the pinned set of type typ with pins
    """
    # assert len(typ) == len(pins)
    return [f for f in flags if projection(f, typ) == pins]


def pincode(poset, xind, zind):
    """create the parity matrices for X and Z checks
    defined by the poset and with xind and zind length
    pinned sets
    """
    poset.makecomplete()
    # print(poset.transitions)
    assert (xind + zind) < poset.length
    flags = poset.get_flags()
    # print(flags)
    xpsets = sum([ps for ps in poset.get_all_pinned_sets(xind).values()], [])
    zpsets = sum([ps for ps in poset.get_all_pinned_sets(zind).values()], [])
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


def reduced_chain_complex(poset, typ):
    """build the reduced chain complex according to type typ
    assumes the code is x=|typ| z = D-x.
    """
    assert poset.iscomplete
    typ_length = len(typ)
    compl_typ = tuple(set(range(poset.length)) - set(typ))
    compl_length = len(compl_typ)
    reduced_qubits = poset.get_all_pinned_sets_with_bound(compl_length)[compl_typ]
    typ_checks_ps = poset.get_all_pinned_sets(typ_length)[typ]
    co_typ_length = poset.length - 1 - typ_length
    co_checks_ps = []
    for co_typ in combinations(compl_typ, co_typ_length):
        co_checks_ps += poset.get_all_pinned_sets(co_typ_length)[co_typ]
    nqubits = len(reduced_qubits)
    nxchecks = len(typ_checks_ps)
    nzchecks = len(co_checks_ps)
    xmat = np.zeros((nxchecks, nqubits), dtype='uint8')
    zmat = np.zeros((nzchecks, nqubits), dtype='uint8')
    for j, reduced_q in enumerate(reduced_qubits):
        for k, co_check in enumerate(co_checks_ps):
            if all(flag in co_check for flag in reduced_q):
                zmat[k, j] = 1
        for k, typ_check in enumerate(typ_checks_ps):
            if any(flag in typ_check for flag in reduced_q):
                xmat[k, j] = 1
    return xmat, zmat, reduced_qubits


def reduced_vec_to_flag_vec(poset, reduced_vec, reduced_qubits):
    """construct a vector for the pin code
    from a vector in the reduced chain complex.
    """
    flags = poset.get_flags()
    nqubits = len(flags)
    flag_vec = np.zeros((nqubits,), dtype='uint8')
    for j in np.nonzero(reduced_vec)[0]:
        for flag in reduced_qubits[j]:
            flag_vec[flags.index(flag)] = 1
    return flag_vec
