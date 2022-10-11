#!/usr/bin/env python3

import math
import numpy as np
import itertools as it
from icecream import ic
from scipy.special import binom

import combischeme_output


def get_level_of_index(index, lmax):
    if index == 0:
        return 0
# cf. https://python-programs.com/python-program-to-find-position-of-rightmost-set-bit/
# def getFirstSetBitPosition(numb):
    # Calculate and the value of log2(n&-n)+1 which gives the first set bit position
    # of the given number and store it in a variable say result_pos.
    result_pos = math.log2(index & -index)
    # Return the value of result_pos(Which is the position of the first set bit).
    return lmax - int(result_pos)


def partition_integer(integer_to_partition):
    """compute partition of positive integer
    -- cf. accel_ast at https://jeromekelleher.net/generating-integer-partitions.html"""
    a = [0 for i in range(integer_to_partition + 1)]
    k = 1
    y = integer_to_partition - 1
    while k != 0:
        x = a[k - 1] + 1
        k -= 1
        while 2 * x <= y:
            a[k] = x
            y -= x
            k += 1
        l = k + 1
        while x <= y:
            a[k] = x
            a[l] = y
            yield a[:k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        yield a[:k + 1]
    return a


def partition_integer_in_num_partitions(integer_to_partition, num_partitions):
    assert (num_partitions <= integer_to_partition)
    partition_list = list(partition_integer(integer_to_partition))
    partition_list = [p for p in partition_list if len(p) == num_partitions]
    return partition_list


def partition_integer_in_num_partitions_with_zeros(integer_to_partition, num_partitions):
    partition_list = []
    shorter = partition_integer(
        integer_to_partition)
    for partition in shorter:
        if len(list(partition)) > num_partitions:
            continue
        # add zeros to partition to make it of length num_partitions
        partition = list(partition)
        partition.extend([0]*(num_partitions-len(partition)))
        # ic(partition)
        permutations = [list(p) for p in it.permutations(partition)]
        partition_list.extend(permutations)
    return partition_list


def get_downward_closed_set_from_level_vector(level_vector):
    subs = [list(range(0, x + 1)) for x in level_vector]
    down_set = it.product(*subs)
    return down_set


def get_min_level_sum(lmin, lmax):
    maxInd = np.argmax(np.array(lmax)-np.array(lmin))
    lm = np.array(lmin)
    lm[maxInd] = lmax[maxInd]
    return lm.sum()


def get_num_dof_of_full_grid(level_vector, boundary):
    for b in boundary:
        assert (b == 2)
    return np.prod([2**l + 1 for l in level_vector])


def get_num_dof_of_subspace(level_vector, boundary):
    for b in boundary:
        assert (b == 2)
    return np.prod([2**(l-1) if l > 0 else 2 for l in level_vector])

# computes the active set by minimum level difference
# todo also implement more sensible schemes like the tilted plane by Christoph Kowitz


def compute_active_set(lmin, lmax):
    listOfRanges = [list(range(0, lmax[i]+1))
                    for i in range(len(lmax))]
    listOfAllGridsUpToLmax = list(it.product(*listOfRanges))
    diagonalIndex = 0
    levelSum = get_min_level_sum(lmin, lmax)+diagonalIndex
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
        for q in range(dim):
            coeff = (-1)**q * binom(dim-1, q)
            levelSum = sum(lmin) + firstLevelDifference - q
            # ic(coeff, levelSum)
            for offset in partition_integer_in_num_partitions_with_zeros(firstLevelDifference-q, dim):
                # ic(offset)
                grid = np.array(lmin) + np.array(offset)
                if (np.sum(grid) == levelSum and (np.array(grid) >= np.array(lmin)).all()):
                    combination_dictionary[tuple(grid)] = coeff

    else:
        # algorithm that can be derived by Hamming distance
        # (cf. Harding 2016, https://link.springer.com/content/pdf/10.1007%2F978-3-319-28262-6_4.pdf)
        ic("non-binomial")
        ic("warning: no other active set implemented yet, will give same result as binomial but take longer!")
        for l in compute_active_set(lmin, lmax):
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


def get_combination_dictionary_from_file(filename):
    combination_dictionary = {}
    data = combischeme_output.read_data_from_json(filename)
    for grid in data:
        combination_dictionary[tuple(grid['level'])] = grid['coeff']
    return combination_dictionary


class CombinationScheme():
    def get_lmax(self):
        return self._lmax

    def get_lmin(self):
        return self._lmin

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

    def get_total_num_points_combi(self):
        total_num_points = 0
        for level in self.get_levels_of_nonzero_coefficient():
            total_num_points += get_num_dof_of_full_grid(
                level, self._boundary_points)
        return total_num_points

    def get_num_grids_per_level_sum(self):
        num_grids_per_level_sum = {}
        for level in self.get_levels_of_nonzero_coefficient():
            levelSum = np.sum([l for l in level])
            if levelSum in num_grids_per_level_sum:
                num_grids_per_level_sum[levelSum] += 1
            else:
                num_grids_per_level_sum[levelSum] = 1
        return num_grids_per_level_sum

    def get_combination_dictionary(self):
        return self._combination_dictionary

    def get_levels_of_nonzero_coefficient(self):
        return self._combination_dictionary.keys()

    def get_nonzero_coefficients(self):
        return self._combination_dictionary.values()


class CombinationSchemeFromMaxLevel(CombinationScheme):
    def __init__(self, lmax, lmin=None, boundary_points=None):
        self._lmax = lmax
        if lmin == None:
            self._lmin = [1]*len(self.lmax)
        else:
            self._lmin = lmin
        if boundary_points == None:
            self._boundary_points = [2]*len(self._lmax)
        else:
            self._boundary_points = boundary_points
            raise NotImplemented
        assert (len(self._lmin) == len(self._lmax))
        for i in range(len(lmax)):
            assert (lmin[i] <= lmax[i])
        self._combination_dictionary = compute_combination_dictionary(
            lmin, lmax)
        assert (self._combination_dictionary is not None)


class CombinationSchemeFromCombinationDictionary(CombinationScheme):
    def __init__(self, dictionary, boundary_points=None):
        self._combination_dictionary = dictionary
        assert (self._combination_dictionary is not None)
        self._lmax = np.max(
            list(self.get_levels_of_nonzero_coefficient()), axis=0)
        self._lmin = np.min(
            list(self.get_levels_of_nonzero_coefficient()), axis=0)
        ic(self._lmax, self._lmin)
        assert (len(self._lmin) == len(self._lmax))
        for i in range(len(self._lmax)):
            assert (self._lmin[i] <= self._lmax[i])
        if boundary_points == None:
            self._boundary_points = [2]*len(self._lmax)
        else:
            self._boundary_points = boundary_points
            raise NotImplemented


class CombinationSchemeFromFile(CombinationSchemeFromCombinationDictionary):
    def __init__(self, filename: str, boundary_points=None):
        dictionary = get_combination_dictionary_from_file(
            filename)
        super().__init__(dictionary, boundary_points)


def write_scheme_to_json(scheme: CombinationScheme):
    # assuming double data type = 8 bytes = 64 bit
    mem = (scheme.get_total_num_points_combi()*8)
    # ic(readable_bytes(mem))
    filename = 'scheme_' + combischeme_output.readable_bytes(mem) + '.json'
    combischeme_output.write_scheme_dictionary_to_json(
        scheme.get_combination_dictionary(), filename)
