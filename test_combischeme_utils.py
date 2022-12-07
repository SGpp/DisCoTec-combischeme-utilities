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
    assert (len(scheme.get_nonzero_coefficients()) > 1)


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


def test_necessary_sparse_grid_spaces(dim=3):
    lmin = [1]*dim
    lmax = [6]*dim
    scheme = combischeme_utils.CombinationSchemeFromMaxLevel(lmax, lmin)
    subspaces = scheme.get_sparse_grid_spaces()
    # this should be a downward closed set up to a reduced level sum
    necessary_subspaces = scheme.get_necessary_sparse_grid_spaces()
    assert len(necessary_subspaces) < len(subspaces)
    assert (tuple(lmin) in list(necessary_subspaces))
    assert (tuple([0]*len(lmin)) in list(necessary_subspaces))
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
    if dim == 3:
        assert num_sg_dof == 1505
    if dim == 6:
        assert num_sg_dof == 159489

    other_scheme = combischeme_utils.CombinationSchemeFromCombinationDictionary(
        scheme.get_combination_dictionary())
    other_subspaces = scheme.get_sparse_grid_spaces()
    other_necessary_subspaces = other_scheme.get_necessary_sparse_grid_spaces()
    assert other_subspaces == subspaces
    assert other_necessary_subspaces == necessary_subspaces


def test_get_num_dof():
    boundary = [2]*3
    level_vectors = [[2, 2, 2], [2, 2, 3], [2, 3, 2],  [3, 2, 2]]
    # fg dof should be 5^3 + 3 * 5^2*9 = 800
    num_fg_dof = combischeme_utils.get_num_dof_of_full_grids(
        level_vectors, boundary)
    assert num_fg_dof == 800
    # sg dof should be 2^3 + 3*2^2*4 = 56
    num_sg_dof = combischeme_utils.get_num_dof_of_subspaces(
        level_vectors, boundary)
    assert num_sg_dof == 56


def test_split_scheme():
    dim = 4
    lmin = [1]*dim
    lmax = [6]*dim
    boundary = [1]*dim
    scheme = combischeme_utils.CombinationSchemeFromMaxLevel(lmax, lmin)
    scheme1, scheme2 = combischeme_utils.split_scheme_by_level_sum(
        scheme)

    # check if the split is correct
    for level_vector in scheme.get_levels_of_nonzero_coefficient():
        found = 0
        if level_vector in scheme1.get_levels_of_nonzero_coefficient():
            found += 1
        if level_vector in scheme2.get_levels_of_nonzero_coefficient():
            found += 1
        assert found == 1

    # compute necessary sg dofs before
    sg_spaces_initial = scheme.get_necessary_sparse_grid_spaces()
    ic(combischeme_utils.get_num_dof_of_subspaces(
        sg_spaces_initial, boundary))
    num_dof_sg_initial = combischeme_utils.get_num_dof_of_subspaces(
        sg_spaces_initial, boundary)
    assert (num_dof_sg_initial == 3072)

    # compute necessary sg dofs after
    subspaces1 = scheme1.get_necessary_sparse_grid_spaces()
    num_dof_sg_1 = combischeme_utils.get_num_dof_of_subspaces(
        subspaces1, boundary)
    subspaces2 = scheme2.get_necessary_sparse_grid_spaces()
    num_dof_sg_2 = combischeme_utils.get_num_dof_of_subspaces(
        subspaces2, boundary)
    ic(num_dof_sg_1, num_dof_sg_2)
    assert (num_dof_sg_1 == 2048)
    assert (num_dof_sg_2 == 2048)
    assert subspaces1.union(subspaces2) == sg_spaces_initial

    # compute conjoint sg dofs
    conjoint_reduced_subspaces = subspaces1.intersection(subspaces2)
    num_conjoint_reduced_dof = combischeme_utils.get_num_dof_of_subspaces(
        conjoint_reduced_subspaces, boundary)
    ic(num_conjoint_reduced_dof)
    assert (num_conjoint_reduced_dof == 1024)
    assert num_conjoint_reduced_dof < num_dof_sg_1
    assert num_conjoint_reduced_dof < num_dof_sg_2

    # compare to full subspaces
    all_subspaces = scheme.get_sparse_grid_spaces()
    all_subspaces1 = scheme1.get_sparse_grid_spaces()
    all_subspaces2 = scheme2.get_sparse_grid_spaces()
    assert all_subspaces1.union(all_subspaces2) == all_subspaces
    conjoint_all_subspaces = all_subspaces1.intersection(all_subspaces2)
    assert conjoint_all_subspaces == conjoint_reduced_subspaces


def test_split_scheme_metis():
    dim = 4
    lmin = [1]*dim
    lmax = [6]*dim
    boundary = [1]*dim
    scheme = combischeme_utils.CombinationSchemeFromMaxLevel(lmax, lmin)

    # TODO it would be nicer if this was installed in the CI
    try:
        import metis
    except (ImportError, RuntimeError):
        print("Skipping test_split_scheme_metis because metis is not installed")
        return

    # use metis split
    schemes_metis = scheme.split_scheme_metis(2)
    assert (len(schemes_metis) == 2)
    scheme1_metis = schemes_metis[0]
    scheme2_metis = schemes_metis[1]

    # TODO check that all are assigned

    # compute necessary sg dofs before
    sg_spaces_initial = scheme.get_necessary_sparse_grid_spaces()
    ic(combischeme_utils.get_num_dof_of_subspaces(
        sg_spaces_initial, boundary))

    # compute necessary sg dofs after
    subspaces1 = scheme1_metis.get_necessary_sparse_grid_spaces()
    num_dof_sg_1 = combischeme_utils.get_num_dof_of_subspaces(
        subspaces1, boundary)
    subspaces2 = scheme2_metis.get_necessary_sparse_grid_spaces()
    num_dof_sg_2 = combischeme_utils.get_num_dof_of_subspaces(
        subspaces2, boundary)
    ic(num_dof_sg_1, num_dof_sg_2)
    assert subspaces1.union(subspaces2) == sg_spaces_initial

    # check if they are the same number as for level-sum split
    assert (num_dof_sg_1 == 2048)
    assert (num_dof_sg_2 == 2048)

    # TODO check conjoint spaces /dofs


def test_integration_create_split_assign(dim=4):
    lmin = [1]*dim
    lmax = [6]*dim
    scheme = combischeme_utils.CombinationSchemeFromMaxLevel(lmax, lmin)
    scheme1, scheme2 = combischeme_utils.split_scheme_by_level_sum(scheme)
    num_process_groups1 = 16
    num_process_groups2 = 13
    assignment1, assigned_FG_size1 = combischeme_utils.assign_combischeme_to_groups(
        scheme1, num_process_groups1)
    assignment2, assigned_FG_size2 = combischeme_utils.assign_combischeme_to_groups(
        scheme2, num_process_groups2)
    assert (len(assignment1) == num_process_groups1)
    assert (len(assignment2) == num_process_groups2)
    assert np.max(assigned_FG_size1) < np.min(assigned_FG_size1)*1.2
    assert np.max(assigned_FG_size1) < np.min(assigned_FG_size2)

    combischeme_output.write_assignment_to_json(
        assignment1, "scheme_test_split_1_"+str(num_process_groups1)+"groups.json")

    combischeme_output.write_assignment_to_json(
        assignment2, "scheme_test_split_2_"+str(num_process_groups2)+"groups.json")
