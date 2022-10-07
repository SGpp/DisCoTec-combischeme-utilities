#!/usr/bin/env python3

import argparse
import numpy as np
from icecream import ic
import json
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
        default=[4,5,6],
    )
    args = parser.parse_args()
    # access CLI options
    lmin = args.lmin
    lmax = args.lmax

    ic(lmin,lmax)

    scheme = combischeme_utils.CombinationScheme(lmax, lmin)

    ic(scheme.get_num_component_grids())
    ic(scheme.get_num_grids_per_level_sum())
    ic(scheme.get_total_num_points_combi())

    # minimum memory requirement of full grids in scheme in bytes
    mem = (scheme.get_total_num_points_combi()*8)
    ic(combischeme_output.readable_bytes(mem))

    schemeList = []
    for key, value in scheme.get_combination_dictionary().items():
        # ic(assignment[group_no])
        schemeList += [{"coeff": value, "level": list(key)}]

    # ic(schemeList)
    jsonString = json.dumps(schemeList)  # , indent=0)

    with open('scheme_' + combischeme_output.readable_bytes(mem) + '.json', 'w') as f:
        f.write(jsonString)
