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
    ic(totalNumPointsCombi, totalNumPointsCombi/1e13)

    levels = list(scheme.get_levels_of_nonzero_coefficient())
    levels.sort(key=lambda x: combischeme_utils.get_num_dof_of_full_grid(
        x, [2]*dim), reverse=True)
    ic(levels[:5])

    assignment = []
    for i in range(num_process_groups):
        assignment.append({})
    assigned_FG_size = [0.]*num_process_groups
    # ic(assignment,assigned_FG_size)
    nextIndex = 0
    for level in levels:
        assignment[nextIndex][level] = scheme.get_coefficient(level)
        assigned_FG_size[nextIndex] += combischeme_utils.get_num_dof_of_full_grid(level, [
                                                                                  2]*dim)
        # this is where load balancing happens!
        nextIndex = np.argmin(assigned_FG_size)
    ic(assigned_FG_size)
    ic([len(a) for a in assignment])
    assert (sum(assigned_FG_size) == totalNumPointsCombi)
    assert (sum([len(a) for a in assignment]) == len(levels))

    for a in assigned_FG_size:
        assert (a*2*8/1e9/(5**6) < 1.)

    combischeme_output.write_assignment_to_json(
        assignment, filename[:-5]+"_"+str(num_process_groups)+"groups.json")
