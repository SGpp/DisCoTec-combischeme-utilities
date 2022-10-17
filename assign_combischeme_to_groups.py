#!/usr/bin/env python3

import argparse
import numpy as np
from icecream import ic
import combischeme_utils
import combischeme_output


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--file_name",
        type=str,
        default="scheme.json",
    )
    parser.add_argument(
        "--num_groups",
        type=int,
        default="64",
    )

    args = parser.parse_args()
    filename = args.file_name
    num_process_groups = args.num_groups
    ic(filename, num_process_groups)

    scheme = combischeme_utils.CombinationSchemeFromFile(filename)

    dim = scheme.get_dimensionality()
    lmax = scheme.get_lmax()
    ic(dim, lmax)

    totalNumPointsCombi = scheme.get_total_num_points_combi()

    ic(scheme.get_num_grids_per_level_sum())
    ic(totalNumPointsCombi, totalNumPointsCombi/1e13,
       combischeme_output.readable_bytes(totalNumPointsCombi*8))

    assignment, assigned_FG_size = combischeme_utils.assign_combischeme_to_groups(
        scheme, num_process_groups)

    for a in assigned_FG_size:
        assert (a*2*8/1e9/(5**6) < 1.)

    combischeme_output.write_assignment_to_json(
        assignment, filename[:-5]+"_"+str(num_process_groups)+"groups.json")
