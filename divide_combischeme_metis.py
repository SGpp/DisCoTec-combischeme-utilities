#!/usr/bin/env python3

import argparse
import numpy as np
from icecream import ic
import combischeme_output
import combischeme_utils

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--lmin",
        nargs="*",
        type=int,
        default=[1, 1, 1],
    )
    parser.add_argument(
        "--lmax",
        nargs="*",
        type=int,
        default=[4, 5, 6],
    )

    args = parser.parse_args()
    lmin = args.lmin
    lmax = args.lmax

    num_partitions = 2

    scheme = combischeme_utils.CombinationSchemeFromMaxLevel(
        lmax, lmin, boundary_points=[1]*len(lmin))
    boundary = scheme.get_boundary_points()

    # compute sg dofs before
    sg_dof_initial = combischeme_utils.get_num_dof_of_subspaces(
        scheme.get_sparse_grid_spaces(), boundary)
    ic(sg_dof_initial, combischeme_output.readable_bytes(sg_dof_initial*8))

    partition_schemes = scheme.split_scheme_metis(num_partitions)
    assert (len(partition_schemes) == num_partitions)
    scheme1 = partition_schemes[0]
    scheme2 = partition_schemes[1]

    # ic(scheme1.get_combination_dictionary(), scheme2.get_combination_dictionary())

    # # compute sg dofs after
    # subspaces1 = scheme1.get_sparse_grid_spaces()
    # subspaces2 = scheme2.get_sparse_grid_spaces()
    # sg_dof1 = combischeme_utils.get_num_dof_of_subspaces(
    #     subspaces1, boundary)
    # sg_dof2 = combischeme_utils.get_num_dof_of_subspaces(
    #     subspaces2, boundary)
    # ic(sg_dof1, sg_dof2)
    # ic(combischeme_output.readable_bytes(sg_dof1*8),
    #    combischeme_output.readable_bytes(sg_dof2*8))

    from itertools import product

    all_metis_params = {
        "objtype": ['cut', 'vol'],
        "ctype": ['rm', 'shem'],
        "iptype": ['grow', 'random', 'edge', 'node'],
        "rtype": ['fm', 'greedy', 'sep2sided', 'sep1sided'],
        "ncuts": [1, 10, 100, 1000, 10000, 100000],
        "niter": [1, 10, 100, 1000, 10000, 100000, 1000000, 10000000],
        "ufactor": [1, 10, 100, 200],
        "minconn": [True, False]
    }
    # iterate all parameter combinations within a (sensible?) range
    param_names = list(all_metis_params.keys())
    # zip with parameter names in order to get original property
    param_values = list(dict(zip(param_names, x))
                        for x in product(*all_metis_params.values()))
    ic(param_values[0])

    current_best = sg_dof_initial

    for metis_args in param_values:
        partition_schemes = scheme.split_scheme_metis(
            num_partitions, **metis_args)

        # compute conjoint sg dofs
        conjoint_subspaces = combischeme_utils.get_conjoint_subspaces(
            partition_schemes)
        sg_dof_conjoint = combischeme_utils.get_num_dof_of_subspaces(
            conjoint_subspaces, boundary)
        if sg_dof_conjoint < current_best:
            current_best = sg_dof_conjoint
            ic("new best!", metis_args, sg_dof_conjoint,
               combischeme_output.readable_bytes(sg_dof_conjoint*8))
            print(metis_args, sg_dof_conjoint,
                  combischeme_output.readable_bytes(sg_dof_conjoint*8))
