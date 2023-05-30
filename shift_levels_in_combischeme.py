#!/usr/bin/env python3

import argparse
from icecream import ic
import combischeme_utils
import combischeme_output


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "file_name",
        type=str,
        default="scheme.json",
    )
    parser.add_argument(
        "--level",
        nargs="*",
        type=int
    )

    args = parser.parse_args()
    level = args.level
    filename = args.file_name
    ic(filename, level)

    scheme = combischeme_utils.CombinationSchemeFromFile(
        filename, boundary_points=1)

    dim = scheme.get_dimensionality()
    assert (dim == len(level))
    lmaxBefore = scheme.get_lmax()
    ic(dim, lmaxBefore)
    totalNumPointsCombiBefore = scheme.get_total_num_points_combi()

    ic(scheme.get_num_grids_per_level_sum())
    ic(totalNumPointsCombiBefore, totalNumPointsCombiBefore/1e13,
       combischeme_output.readable_bytes(totalNumPointsCombiBefore*8))

    # change every key in combination dictionary by level
    combidictBefore = scheme.get_combination_dictionary()
    data = combischeme_output.read_data_from_json(filename)
    levels = []
    coeffs = []
    groupNos = []
    try:
        for grid in data:
            levels.append(tuple(grid['level']))
            coeffs.append(grid['coeff'])
            groupNos.append(grid['group_no'])
        assert (len(levels) == len(coeffs) == len(groupNos))
        assignmentBefore = [{}.copy() for _ in range(max(groupNos) + 1)]
        for i in range(len(levels)):
            assignmentBefore[groupNos[i]][levels[i]] = coeffs[i]
        assignmentAfter = [{}.copy() for _ in range(max(groupNos) + 1)]
        for groupNo in range(len(assignmentBefore)):
            for key in assignmentBefore[groupNo]:
                assignmentAfter[groupNo][tuple([l + level[i]
                                                for i, l in enumerate(key)])] = combidictBefore[key]
        ic(len(levels), max([len(a) for a in assignmentAfter]), min(
            [len(a) for a in assignmentAfter]))
    except KeyError as k:
        print("no group assignments?")

    combidictAfter = {}
    for key in combidictBefore:
        combidictAfter[tuple([l + level[i]
                             for i, l in enumerate(key)])] = combidictBefore[key]

    schemeAfter = combischeme_utils.CombinationSchemeFromCombinationDictionary(
        combidictAfter, boundary_points=scheme.get_boundary_points())
    totalNumPointsCombiAfter = schemeAfter.get_total_num_points_combi()
    ic(totalNumPointsCombiAfter, totalNumPointsCombiAfter/1e13,
       combischeme_output.readable_bytes(totalNumPointsCombiAfter*8))
    assert (totalNumPointsCombiAfter ==
            totalNumPointsCombiBefore * 2**(sum(level)))
    lmaxAfter = schemeAfter.get_lmax()
    ic(lmaxAfter)
    assert (all(lmaxAfter == [l + level[i] for i, l in enumerate(lmaxBefore)]))

    if len(groupNos) > 0:
        combischeme_output.write_assignment_to_json(assignmentAfter, "scheme_large_" + '-'.join(
            [str(l) for l in lmaxAfter]) + "_" + format(len(assignmentAfter), '05d') + "groups.json")
    else:
        combischeme_utils.write_scheme_to_json(
            schemeAfter, "scheme_large_" + '-'.join([str(l) for l in lmaxAfter]) + ".json")
