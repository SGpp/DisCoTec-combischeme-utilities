#!/usr/bin/env python3

import numpy as np
from icecream import ic
import json
import combischeme_utils
import combischeme_output

#TODO this could be cl-input
lmin = [2]*6
lmax = [4]*6

scheme = combischeme_utils.CombinationScheme(lmax, lmin)

numGridsOfSize = {}

totalNumPointsCombi = 0
for key in scheme.get_levels_of_nonzero_coefficient():
    totalNumPointsCombi += np.prod([2**l + 1 for l in key])
    levelSum=np.sum([l for l in key])
    if levelSum in numGridsOfSize:
        numGridsOfSize[levelSum] += 1
    else:
        numGridsOfSize[levelSum] = 1

ic(scheme.get_num_component_grids)
ic(numGridsOfSize)
ic(totalNumPointsCombi, totalNumPointsCombi/1e13)

# minimum memory requirement of full grids in scheme in bytes
mem = (totalNumPointsCombi*8)
ic(combischeme_output.readable_bytes(mem))

schemeList = []
for key, value in scheme.get_combination_dictionary().items():
    # ic(assignment[group_no])
    schemeList += [{"coeff": value, "level": list(key)}]

# ic(schemeList)
jsonString = json.dumps(schemeList)#, indent=0)

with open('scheme_'+ combischeme_output.readable_bytes(mem) + '.json', 'w') as f:
    f.write(jsonString)
