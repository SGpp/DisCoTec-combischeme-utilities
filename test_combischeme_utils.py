#!/usr/bin/env python3
import combischeme_utils
from icecream import ic
import numpy as np
import pytest

# run this test file with `pytest test_combischeme_utils.py``


def validate_combischeme(scheme: combischeme_utils.CombinationScheme):
    assert (scheme.get_levels_of_nonzero_coefficient() is not None)
    assert (len(list(scheme.get_levels_of_nonzero_coefficient())
            [0]) == scheme.get_dimensionality())
    assert (scheme.get_nonzero_coefficients() is not None)
    ic(scheme.get_nonzero_coefficients())
    assert (not (0. in scheme.get_nonzero_coefficients()))
    sum_coefficients = sum(scheme.get_nonzero_coefficients())
    assert (sum_coefficients == 1.)


def test_create_regular_combischeme(dim=3):
    lmin = [1]*dim
    lmax = [6]*dim
    scheme = combischeme_utils.CombinationScheme(lmax, lmin)
    print(np.sum(scheme.get_nonzero_coefficients))
    validate_combischeme(scheme)


def test_create_nonregular_combischeme(dim=3):
    lmin = [1]*dim
    lmax = [6]*dim
    lmax[0] = 5
    scheme = combischeme_utils.CombinationScheme(lmax, lmin)
    validate_combischeme(scheme)
