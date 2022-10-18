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
    parser.add_argument(
        "--num_groups",
        nargs=2,
        type=int,
        default=[32,16],
    )
    args = parser.parse_args()
    # access CLI options
    lmin = args.lmin
    lmax = args.lmax
    num_process_groups = args.num_groups

    ic(lmin, lmax)

    scheme = combischeme_utils.CombinationSchemeFromMaxLevel(lmax, lmin)

    dim = scheme.get_dimensionality()
    lmax = scheme.get_lmax()
    ic(dim, lmax)

    num_grids_per_level_sum = scheme.get_num_grids_per_level_sum()
    highest_level_sum = int(np.max(list(num_grids_per_level_sum.keys())))
    ic(num_grids_per_level_sum, highest_level_sum)

    scheme1, scheme2 = combischeme_utils.split_scheme_by_level_sum(
        scheme)

    assignment1, assigned_FG_size1 = combischeme_utils.assign_combischeme_to_groups(
        scheme1, num_process_groups[0])

    assignment2, assigned_FG_size2 = combischeme_utils.assign_combischeme_to_groups(
        scheme2, num_process_groups[1])

    combischeme_output.write_assignment_to_json(
        assignment1, "scheme_large_" + '-'.join([str(l) for l in lmax]) + "_split1_"+str(num_process_groups[0])+"groups.json")
    combischeme_output.write_assignment_to_json(
        assignment2, "scheme_large_" + '-'.join([str(l) for l in lmax]) + "_split2_"+str(num_process_groups[1])+"groups.json")

    # here goes all the diagnostic output, which takes long and you can abort the script if you're not interested
    ic("diagnostic output")
    # minimum memory requirement of full grids in scheme in bytes
    ic(combischeme_output.readable_bytes(scheme.get_total_num_points_combi()*8))
    # divided memory requirement
    ic(combischeme_output.readable_bytes(scheme1.get_total_num_points_combi()*8))
    ic(combischeme_output.readable_bytes(scheme2.get_total_num_points_combi()*8))

    # ic(scheme.get_total_num_points_combi())
    ic(scheme.get_num_component_grids())
    ic(scheme.get_num_grids_per_level_sum())

    # compute sg dofs before
    sg_dof_initial = combischeme_utils.get_num_dof_of_subspaces(
        scheme.get_necessary_sparse_grid_spaces(), [2]*scheme.get_dimensionality())
    ic(sg_dof_initial, combischeme_output.readable_bytes(sg_dof_initial*8))

    # compute sg dofs after
    subspaces1 = scheme1.get_necessary_sparse_grid_spaces()
    subspaces2 = scheme2.get_necessary_sparse_grid_spaces()
    sg_dof1 = combischeme_utils.get_num_dof_of_subspaces(
        subspaces1, [2]*scheme.get_dimensionality())
    sg_dof2 = combischeme_utils.get_num_dof_of_subspaces(
        subspaces2, [2]*scheme.get_dimensionality())
    ic(sg_dof1, sg_dof2)
    ic(combischeme_output.readable_bytes(sg_dof1*8),
       combischeme_output.readable_bytes(sg_dof2*8))

    # compute conjoint sg dofs
    conjoint_subspaces = subspaces1.intersection(subspaces2)
    sg_dof_conjoint =combischeme_utils.get_num_dof_of_subspaces(
        conjoint_subspaces, [2]*scheme.get_dimensionality())
    ic(sg_dof_conjoint, combischeme_output.readable_bytes(sg_dof_conjoint*8))
