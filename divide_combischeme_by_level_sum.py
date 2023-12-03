#!/usr/bin/env python3

import argparse
import numpy as np
from icecream import ic
import combischeme_output
import combischeme_utils

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--file_name",
        type=str,
        default="scheme.json",
    )

    args = parser.parse_args()
    filename = args.file_name

    scheme = combischeme_utils.CombinationSchemeFromFile(filename)

    dim = scheme.get_dimensionality()
    boundary = [1]*dim
    lmax = scheme.get_lmax()
    ic(dim, lmax)

    num_grids_per_level_sum = scheme.get_num_grids_per_level_sum()
    highest_level_sum = int(np.max(list(num_grids_per_level_sum.keys())))
    ic(num_grids_per_level_sum, highest_level_sum)

    scheme1, scheme2 = combischeme_utils.split_scheme_by_level_sum(scheme)
    combischeme_utils.write_scheme_to_json(
        scheme1, filename[:-5]+"_split1.json")
    combischeme_utils.write_scheme_to_json(
        scheme2, filename[:-5]+"_split2.json")

    # compute sg dofs before
    sg_dof_initial = combischeme_utils.get_num_dof_of_subspaces(
        scheme.get_necessary_sparse_grid_spaces(), boundary)
    ic(sg_dof_initial, combischeme_output.readable_bytes(sg_dof_initial*8))

    # compute sg dofs after
    subspaces1 = scheme1.get_necessary_sparse_grid_spaces()
    subspaces2 = scheme2.get_necessary_sparse_grid_spaces()
    sg_dof1 = combischeme_utils.get_num_dof_of_subspaces(
        subspaces1, boundary)
    sg_dof2 = combischeme_utils.get_num_dof_of_subspaces(
        subspaces2, boundary)
    ic(sg_dof1, sg_dof2)
    ic(combischeme_output.readable_bytes(sg_dof1*8),
       combischeme_output.readable_bytes(sg_dof2*8))

    # compute conjoint sg dofs
    conjoint_subspaces = subspaces1.intersection(subspaces2)
    sg_dof_conjoint = combischeme_utils.get_num_dof_of_subspaces(
        conjoint_subspaces, boundary)
    ic(sg_dof_conjoint, combischeme_output.readable_bytes(sg_dof_conjoint*8))
