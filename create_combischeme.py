#!/usr/bin/env python3

import argparse
from icecream import ic
import combischeme_utils
import combischeme_output


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
    # access CLI options
    lmin = args.lmin
    lmax = args.lmax

    ic(lmin, lmax)

    scheme = combischeme_utils.CombinationSchemeFromMaxLevel(
        lmax, lmin, boundary_points=[1]*len(lmin))

    combischeme_utils.write_scheme_to_json(
        scheme)

    # minimum memory requirement of full grids in scheme in bytes
    mem = (scheme.get_total_num_points_combi()*8)
    ic(combischeme_output.readable_bytes(mem))

    ic(scheme.get_num_component_grids())
    ic(scheme.get_num_grids_per_level_sum())
    ic(scheme.get_total_num_points_combi())
    # ic(scheme.get_total_num_points_sparse())

    # compute sg dofs
    sg_dof_initial = combischeme_utils.get_num_dof_of_subspaces(
        scheme.get_necessary_sparse_grid_spaces(), scheme.get_boundary_points())
    ic(sg_dof_initial, combischeme_output.readable_bytes(sg_dof_initial*8))
