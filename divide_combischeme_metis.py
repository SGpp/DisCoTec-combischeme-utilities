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

    # compute minimum memory requirement of full grids in scheme in bytes
    ic(scheme.get_total_num_points_combi(), combischeme_output.readable_bytes(
        scheme.get_total_num_points_combi()*8))
    # compute sg dofs before
    sg_dof_initial = combischeme_utils.get_num_dof_of_subspaces(
        scheme.get_sparse_grid_spaces(), boundary)
    ic(sg_dof_initial, combischeme_output.readable_bytes(sg_dof_initial*8))

    target_partition_weights = [1./num_partitions]*num_partitions
    partition_schemes = scheme.split_scheme_metis(
        num_partitions, tpwgts=target_partition_weights)
    assert (len(partition_schemes) == num_partitions)

    for i in range(num_partitions):
        combischeme_utils.write_scheme_to_json(
            partition_schemes[i], "scheme_large_" + '-'.join([str(l) for l in lmax]) + "_part" + str(i) + ".json")

    scheme1 = partition_schemes[0]
    scheme2 = partition_schemes[1]

    # ic(scheme1.get_combination_dictionary(), scheme2.get_combination_dictionary())

    ic("diagnostic output")
    # divided memory requirement
    ic(scheme1.get_total_num_points_combi(), combischeme_output.readable_bytes(
        scheme1.get_total_num_points_combi()*8))
    ic(scheme2.get_total_num_points_combi(), combischeme_output.readable_bytes(
        scheme2.get_total_num_points_combi()*8))

    # compute sg dofs after
    subspaces1 = scheme1.get_sparse_grid_spaces()
    subspaces2 = scheme2.get_sparse_grid_spaces()
    sg_dof1 = combischeme_utils.get_num_dof_of_subspaces(
        subspaces1, boundary)
    sg_dof2 = combischeme_utils.get_num_dof_of_subspaces(
        subspaces2, boundary)
    ic(sg_dof1, sg_dof2)
    ic(combischeme_output.readable_bytes(sg_dof1*8),
       combischeme_output.readable_bytes(sg_dof2*8))

    # compute conjoint sg dofs
    conjoint_subspaces = combischeme_utils.get_conjoint_subspaces(
        partition_schemes)
    sg_dof_conjoint = combischeme_utils.get_num_dof_of_subspaces(
        conjoint_subspaces, boundary)
    ic(sg_dof_conjoint)
    ic(combischeme_output.readable_bytes(sg_dof_conjoint*8))

    # output parallelization info
    max_process_group_size = 2**sum(lmin)
    sg_memory_per_process = max(sg_dof1, sg_dof2) / \
        max_process_group_size * 8 / 1e9
    extra_sg_memory_per_process = sg_dof_conjoint / max_process_group_size * 8 / 1e9
    available_memory_per_process = 2.
    memory_left_per_process = available_memory_per_process - \
        sg_memory_per_process - extra_sg_memory_per_process
    fg_memory_to_distribute = scheme.get_total_num_points_combi() / \
        max_process_group_size * 8 / 1e9
    print("if running on the (maximum) process group size of " + str(max_process_group_size) +
          ", this scenario will use " + str(sg_memory_per_process) +
          " GB per process for the sparse grid data; and " + str(extra_sg_memory_per_process) +
          " for the extra sparse grid data. Under optimal conditions, one would need " +
          str(fg_memory_to_distribute/memory_left_per_process) + " process groups if " +
          str(available_memory_per_process) + " GB main memory would be available per core.")
