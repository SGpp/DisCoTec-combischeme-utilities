#!/usr/bin/env python3
import combischeme_utils
import combischeme_output
from icecream import ic
import numpy as np
import scipy
import pytest
import time

# run this test file with `pytest test_combischeme_utils.py``


def validate_combischeme(scheme: combischeme_utils.CombinationScheme):
    assert (scheme.get_levels_of_nonzero_coefficient() is not None)
    assert (len(list(scheme.get_levels_of_nonzero_coefficient())
            [0]) == scheme.get_dimensionality())
    ic(scheme.get_combination_dictionary())
    assert (scheme.get_nonzero_coefficients() is not None)
    ic(scheme.get_nonzero_coefficients())
    assert (not (0. in scheme.get_nonzero_coefficients()))
    sum_coefficients = sum(scheme.get_nonzero_coefficients())
    assert (sum_coefficients == 1.)


def test_create_regular_combischeme(dim=3):
    lmin = [1]*dim
    lmax = [6]*dim
    scheme = combischeme_utils.CombinationSchemeFromMaxLevel(lmax, lmin)
    print(np.sum(scheme.get_nonzero_coefficients))
    for d in range(dim):
        corner = lmin.copy()
        corner[d] = lmax[d]
        assert (tuple(corner) in list(
            scheme.get_levels_of_nonzero_coefficient()))
        for d2 in range(dim):
            corner_neighbor = corner.copy()
            corner_neighbor[d2] += 1
            assert (tuple(corner_neighbor) not in list(
                scheme.get_levels_of_nonzero_coefficient()))
    validate_combischeme(scheme)


def test_create_non_regular_combischeme(dim=3):
    lmin = [1]*dim
    lmax = [6]*dim
    lmax[0] = 5
    scheme = combischeme_utils.CombinationSchemeFromMaxLevel(lmax, lmin)
    validate_combischeme(scheme)


def test_create_degenerate_combischeme(dim=3):
    lmin = [1]*dim
    lmax = [6]*dim
    lmax[0] = lmin[0]
    scheme = combischeme_utils.CombinationSchemeFromMaxLevel(lmax, lmin)
    validate_combischeme(scheme)
    assert(len(scheme.get_nonzero_coefficients()) > 1)


def test_create_combischeme_from_dict():
    dictionary = {(2, 2, 2): 1,
                  (2, 2, 3): -2,
                  (2, 2, 4): 1,
                  (2, 3, 2): -2,
                  (2, 3, 3): 1,
                  (2, 4, 2): 1,
                  (3, 2, 2): -2,
                  (3, 2, 3): 1,
                  (3, 3, 2): 1,
                  (4, 2, 2): 1}
    scheme = combischeme_utils.CombinationSchemeFromCombinationDictionary(
        dictionary)
    validate_combischeme(scheme)


def test_write_combischeme(dim=4):
    lmin = [1]*dim
    lmax = [6]*dim
    scheme = combischeme_utils.CombinationSchemeFromMaxLevel(lmax, lmin)
    combischeme_utils.write_scheme_to_json(scheme)


def test_partition_integer(integer=6):
    partitions = combischeme_utils.partition_integer(integer)
    partitions = list(partitions)
    ic(partitions)
    assert (len(list(partitions[0])) == integer)
    assert ([1]*integer in partitions)
    assert ([integer] in partitions)
    single_partition = combischeme_utils.partition_integer_in_num_partitions(
        integer, 1)
    assert (list(single_partition) == [[integer]])
    single_partition_of_ones = combischeme_utils.partition_integer_in_num_partitions(
        integer, integer)
    assert (list(single_partition_of_ones) == [[1]*integer])
    single_partition_of_ones = combischeme_utils.partition_integer_in_num_partitions(
        integer, integer)
    assert (list(single_partition_of_ones) == [[1]*integer])
    # assert (len(partitions) == scipy.special.binom(integer, i))


def test_partition_integer_in_num_partitions_with_zeros(integer=6):
    filled_partitions = combischeme_utils.partition_integer_in_num_partitions_with_zeros(
        integer, integer)
    for i in range(integer):
        zeros = [0]*(integer)
        zeros[i] = integer
        assert zeros in filled_partitions
    assert ([1]*integer in filled_partitions)


def test_compute_active_set(dim=6):
    lmin = [2]*dim
    lmax = [8]*dim
    tic = time.perf_counter()
    active_set1 = combischeme_utils.compute_active_set(lmin, lmax)
    toc = time.perf_counter()
    time1 = toc - tic
    print(f"first active set in {time1:0.4f}seconds")
    ic(active_set1)
    assert active_set1.intersection([(2, 2, 2, 2, 2, 8)]) != set()
    assert active_set1.intersection([(2, 2, 2, 2, 8, 2)]) != set()
    assert active_set1.intersection([(2, 2, 2, 8, 2, 2)]) != set()
    assert active_set1.intersection([(2, 2, 8, 2, 2, 2)]) != set()
    assert active_set1.intersection([(2, 8, 2, 2, 2, 2)]) != set()
    assert active_set1.intersection([(8, 2, 2, 2, 2, 2)]) != set()

    lmin[0] += 1
    tic = time.perf_counter()
    active_set2 = combischeme_utils.compute_active_set(lmin, lmax)
    toc = time.perf_counter()
    time2 = toc - tic
    print(f"second active set in {time2:0.4f}seconds")

    ic(active_set2)
    # assert active_set1.intersection(active_set2) == active_set2
    assert len(active_set1) == len(active_set2) + 1
    assert time1 < time2


def test_necessary_sparse_grid_spaces(dim=6):
    lmin = [1]*dim
    lmax = [6]*dim
    scheme = combischeme_utils.CombinationSchemeFromMaxLevel(lmax, lmin)
    # this should be a downward closed set up to a reduced level sum
    necessary_subspaces = scheme.get_necessary_sparse_grid_spaces()
    assert (tuple(lmin) in list(necessary_subspaces))
    for d in range(dim):
        corner = lmin.copy()
        corner[d] = lmax[d] - 1
        assert (tuple(corner) in list(necessary_subspaces))
        for d2 in range(dim):
            corner_neighbor = corner.copy()
            corner_neighbor[d2] += 1
            assert (tuple(corner_neighbor) not in list(necessary_subspaces))

    num_sg_dof = combischeme_utils.get_num_dof_of_subspaces(
        necessary_subspaces, boundary=[2]*dim)
    # regression test
    assert num_sg_dof == 159489


def test_get_num_dof():
    boundary = [2]*3
    level_vectors = [[2, 2, 2], [2, 2, 3], [2, 3, 2],  [3, 2, 2]]
    # fg dof should be 5^3 + 3 * 5^2*9 = 800
    num_fg_dof = combischeme_utils.get_num_dof_of_full_grids(level_vectors, boundary)
    assert num_fg_dof ==800
    # sg dof should be 2^3 + 3*2^2*4 = 56
    num_sg_dof = combischeme_utils.get_num_dof_of_subspaces(level_vectors, boundary)
    assert num_sg_dof == 56

