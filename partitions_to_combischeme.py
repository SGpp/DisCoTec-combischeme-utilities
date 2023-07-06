#!/usr/bin/env python3

import argparse
import numpy as np
from math import isclose
from icecream import ic
import combischeme_output
import combischeme_utils

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "partitioned_files",
        nargs="*",
        type=str,
    )
    parser.add_argument(
        "--out_name",
        nargs=1,
        type=str,
        default=None,
    )

    args = parser.parse_args()

    first_scheme = combischeme_utils.CombinationSchemeFromFile(
        args.partitioned_files[0])

    # divide json-dictionary into assigned partitions
    partition_dicts = []

    for file in args.partitioned_files:
        data = combischeme_output.read_data_from_json(file)
        # collect all different process group numbers
        group_numbers = set()
        try:
            for grid in data:
                group_numbers.add(grid['group_no'])

            if len(group_numbers) > 1:
                for group_no in group_numbers:
                    combination_dictionary = {}
                    for grid in data:
                        if grid['group_no'] == group_no:
                            combination_dictionary[tuple(
                                grid['level'])] = grid['coeff']
                    partition_dicts.append(combination_dictionary)
        except KeyError:
            combination_dictionary = {}
            for grid in data:
                combination_dictionary[tuple(grid['level'])] = grid['coeff']
            partition_dicts.append(combination_dictionary)

    partitions = [combischeme_utils.CombinationSchemeFromCombinationDictionary(
        dict) for dict in partition_dicts]

    assignment = [p.get_combination_dictionary() for p in partitions]

    if args.out_name is None:
        args.out_name = "combined_"+args.partitioned_files[0]
    combischeme_output.write_assignment_to_json(assignment, args.out_name)

    # for testing, with checks?
    # combined_scheme = combischeme_utils.CombinationSchemeFromPartitions(partitions)
