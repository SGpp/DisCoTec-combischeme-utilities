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

    first_scheme = combischeme_utils.CombinationSchemeFromFile(args.partitioned_files[0])

    partitions = [combischeme_utils.CombinationSchemeFromFile(file) for file in args.partitioned_files]

    if args.out_name is None:
        args.out_name = "combined_"+args.partitioned_files[0]

    assignment = [p.get_combination_dictionary() for p in partitions]

    combischeme_output.write_assignment_to_json(assignment, args.out_name)

    # for testing, with checks?
    # combined_scheme = combischeme_utils.CombinationSchemeFromPartitions(partitions)