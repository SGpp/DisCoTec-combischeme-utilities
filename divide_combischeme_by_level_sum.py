#!/usr/bin/env python3

import argparse
import numpy as np
from icecream import ic
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
    lmax = scheme.get_lmax()
    ic(dim, lmax)

    num_grids_per_level_sum = scheme.get_num_grids_per_level_sum()
    highest_level_sum = int(np.max(list(num_grids_per_level_sum.keys())))
    ic(num_grids_per_level_sum, highest_level_sum)


    gridsForSystems = list(scheme.get_levels_of_nonzero_coefficient())
    gridsForSystem1 = []
    gridsForSystem2 = []
    ic(len(gridsForSystem1), len(gridsForSystem2), len(gridsForSystems))

    gridsToIterate = gridsForSystems.copy()
    rr = 0
    # # for different splits -- the sum of weights is assumed to be 1.,
    # # if system 1's weight is 0.5, we have an even split.
    weight_system_1 = 0.5
    dim_to_split_at = int(dim/2)
    for level in gridsToIterate:
        max_dims = np.argwhere(level == np.max(level))
        lower_dims_diff_to_lmax = np.array(
            lmax[:dim_to_split_at])-np.array(level)[:dim_to_split_at]
        norm_lower_dims = - \
            np.linalg.norm(lower_dims_diff_to_lmax, ord=0.99)
        higher_dims_diff_to_lmax = np.array(
            lmax[dim_to_split_at:])-np.array(level[dim_to_split_at:])
        norm_higher_dims = - \
            np.linalg.norm(higher_dims_diff_to_lmax, ord=0.99)

        if weight_system_1 * norm_lower_dims > (1.-weight_system_1) * norm_higher_dims:
            gridsForSystems.remove(level)
            gridsForSystem1.append(level)
            continue
        elif weight_system_1 * norm_lower_dims < (1.-weight_system_1) * norm_higher_dims:
            gridsForSystems.remove(level)
            gridsForSystem2.append(level)
            continue
        else:
            if rr % 2 == 0:
                gridsForSystems.remove(level)
                gridsForSystem1.append(level)
            else:
                gridsForSystems.remove(level)
                gridsForSystem2.append(level)
            rr += 1

    assert(len(gridsForSystems) == 0)
    ic(len(gridsForSystem1), len(gridsForSystem2), len(gridsForSystems))

    # build the new combination dictionaries
    dictionary1 = {}
    for level in gridsForSystem1:
        dictionary1[level] = scheme.get_coefficient(level)
    scheme1 = combischeme_utils.CombinationSchemeFromCombinationDictionary(dictionary1)
    combischeme_utils.write_scheme_to_json(scheme1, filename[:-5]+"_split1.json")
    dictionary2 = {}
    for level in gridsForSystem2:
        dictionary2[level] = scheme.get_coefficient(level)
    scheme2 = combischeme_utils.CombinationSchemeFromCombinationDictionary(dictionary2)
    combischeme_utils.write_scheme_to_json(scheme2, filename[:-5]+"_split2.json")
