#!/usr/bin/env python3

import math
import numpy as np
import itertools as it
from icecream import ic
from scipy.special import binom


def get_level_of_index(index, lmax):
    if index == 0:
        return 0
# cf. https://python-programs.com/python-program-to-find-position-of-rightmost-set-bit/
#def getFirstSetBitPosition(numb):
    # Calculate and the value of log2(n&-n)+1 which gives the first set bit position
    # of the given number and store it in a variable say result_pos.
    result_pos = math.log2(index & -index)
    # Return the value of result_pos(Which is the position of the first set bit).
    return lmax - int(result_pos)

def get_downward_closed_set_from_level_vector(level_vector):
    subs = [list(range(0, x + 1)) for x in level_vector]
    down_set = it.product(*subs)
    return down_set


def get_min_level_sum(lmin, lmax):
    maxInd = np.argmax(np.array(lmax)-np.array(lmin))
    lm = np.array(lmin)
    lm[maxInd] = lmax[maxInd]
    return lm.sum()

# computes the active set by minimum level difference
#todo also implement more sensible schemes like the tilted plane by Christoph Kowitz


def compute_active_set(lmin, lmax):
    listOfRanges = [list(range(0, lmax[i]+1))
                    for i in range(len(lmax))]
    listOfAllGridsUpToLmax = list(it.product(*listOfRanges))
    diagonalIndex = 0
    levelSum = get_min_level_sum(lmin,lmax)+diagonalIndex
    s = set()
    for grid in listOfAllGridsUpToLmax:
        if (np.sum(grid) == levelSum and (np.array(grid) >= np.array(lmin)).all()):
            s.add(grid)

    return s


def compute_combination_dictionary(lmin, lmax):
    dim = len(lmin)
    combination_dictionary = {}
    firstLevelDifference = lmax[0] - lmin[0]
    uniformLevelDifference = [(lmax[i] - lmin[i]) ==
                               firstLevelDifference for i in range(dim)]
    if uniformLevelDifference and firstLevelDifference >= dim:
        # implements the standard formula, cf. "Sparse Grids in a Nutshell"
        # (cf. https://link.springer.com/chapter/10.1007/978-3-642-31703-3_3)
        ic("binomial")
        listOfRanges = [list(range(lmin[i], lmax[i]+1))
                        for i in range(len(lmax))]
        ic(listOfRanges)
        for q in range(dim):
            coeff = (-1)**q * binom(dim-1, q)
            levelSum = sum(lmin) + firstLevelDifference - q
            # ic(coeff, levelSum)
            for grid in it.product(*listOfRanges):
                if (np.sum(grid) == levelSum and (np.array(grid) >= np.array(lmin)).all()):
                    combination_dictionary[grid] = coeff

    else:
        # algorithm that can be derived by Hamming distance
        # (cf. Harding 2016, https://link.springer.com/content/pdf/10.1007%2F978-3-319-28262-6_4.pdf)
        ic("non-binomial")
        for l in compute_active_set(lmin,lmax):
            combination_dictionary[l] = 1

        dict_of_subspaces = {}
        for l in combination_dictionary:
            for subspace in get_downward_closed_set_from_level_vector(l):
                if subspace in dict_of_subspaces:
                    dict_of_subspaces[subspace] += 1
                else:
                    dict_of_subspaces[subspace] = 1

        # # remove subspaces which are too much
        while (set(dict_of_subspaces.values()) != set([1])):

            for subspace in dict_of_subspaces:
                currentCount = dict_of_subspaces[subspace]
                if currentCount != 1:
                    diff = currentCount - 1

                    if subspace in combination_dictionary:

                        combination_dictionary[subspace] -= diff
                        if combination_dictionary[subspace] == 0:
                            del combination_dictionary[subspace]

                    else:
                        combination_dictionary[subspace] = -diff

                    for l in get_downward_closed_set_from_level_vector(subspace):
                        dict_of_subspaces[l] -= diff
    return combination_dictionary


class CombinationScheme():
    def __init__(self, lmax, lmin=None):
        self._lmax = lmax
        if lmin == None:
            self._lmin = [1]*len(self.lmax)
        else:
            self._lmin = lmin
        assert (len(self._lmin) == len(self._lmax))
        for i in range(len(lmax)):
            assert (lmin[i] <= lmax[i])
        self._combination_dictionary = compute_combination_dictionary(
            lmin, lmax)
        assert (self._combination_dictionary is not None)

    def get_all_subspaces(self):
        subspacesSet = set()
        for l in self._combination_dictionary:
            for subspace in get_downward_closed_set_from_level_vector(l):
                subspacesSet.add(subspace)
        return subspacesSet

    def get_dimensionality(self):
        return len(self._lmin)

    def get_num_component_grids(self):
        return len(self._combination_dictionary)

    def get_combination_dictionary(self):
        return self._combination_dictionary

    def get_levels_of_nonzero_coefficient(self):
        return self._combination_dictionary.keys()

    def get_nonzero_coefficients(self):
        return self._combination_dictionary.values()

