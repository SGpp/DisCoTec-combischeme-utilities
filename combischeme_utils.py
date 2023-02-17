#!/usr/bin/env python3

from __future__ import annotations

from functools import reduce
import itertools as it
try:
    from icecream import ic
    ic("test ic")
except:
    def ic(*args):
        pass
import math
import networkx as nx
import numpy as np
from operator import and_, or_
from scipy.special import binom

import combischeme_output


def shadows(level_vector_that_maybe_shadows, level_vector_that_is_maybe_shadowed) -> bool:
    return not any(np.greater(level_vector_that_is_maybe_shadowed, level_vector_that_maybe_shadows))


def get_level_of_index(index, lmax) -> int:
    if index == 0:
        return 0
# cf. https://python-programs.com/python-program-to-find-position-of-rightmost-set-bit/
# def getFirstSetBitPosition(numb):
    # Calculate and the value of log2(n&-n)+1 which gives the first set bit position
    # of the given number and store it in a variable say result_pos.
    result_pos = math.log2(index & -index)
    # Return the value of result_pos(Which is the position of the first set bit).
    return lmax - int(result_pos)


def partition_integer(integer_to_partition: int):
    """compute partition of positive integer
    -- cf. accel_ast at https://jeromekelleher.net/generating-integer-partitions.html"""
    a = [0 for i in range(integer_to_partition + 1)]
    k = 1
    y = integer_to_partition - 1
    while k != 0:
        x = a[k - 1] + 1
        k -= 1
        while 2 * x <= y:
            a[k] = x
            y -= x
            k += 1
        l = k + 1
        while x <= y:
            a[k] = x
            a[l] = y
            yield a[:k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        yield a[:k + 1]
    return a


def partition_integer_in_num_partitions(integer_to_partition: int, num_partitions: int) -> list:
    assert (num_partitions <= integer_to_partition)
    partition_list = list(partition_integer(integer_to_partition))
    partition_list = [p for p in partition_list if len(p) == num_partitions]
    return partition_list


def partition_integer_in_num_partitions_with_zeros(integer_to_partition: int, num_partitions: int) -> list:
    partition_list = []
    shorter = partition_integer(
        integer_to_partition)
    for partition in shorter:
        if len(list(partition)) > num_partitions:
            continue
        # add zeros to partition to make it of length num_partitions
        partition = list(partition)
        partition.extend([0]*(num_partitions-len(partition)))
        # ic(partition)
        permutations = [list(p) for p in it.permutations(partition)]
        partition_list.extend(permutations)
    return partition_list


def get_downward_closed_set_from_level_vector(level_vector) -> set:
    subs = [tuple(range(0, x + 1)) for x in level_vector]
    down_set = it.product(*subs)
    return down_set


def get_downward_closed_set_from_level_vectors(level_vectors) -> set:
    down_set = set()
    for level_vector in level_vectors:
        down_set.update(
            get_downward_closed_set_from_level_vector(level_vector))
    return down_set


def is_downward_closed_set(level_vectors) -> bool:
    down_set = get_downward_closed_set_from_level_vectors(level_vectors)
    return down_set.symmetric_difference(set(level_vectors)) == set()


def get_min_level_sum(lmin, lmax) -> int:
    maxInd = np.argmax(np.array(lmax)-np.array(lmin))
    lm = np.array(lmin)
    lm[maxInd] = lmax[maxInd]
    return lm.sum()


def get_num_dof_of_full_grid(level_vector, boundary) -> int:
    return np.prod([2**level_vector[i] - 1 + boundary[i] for i in range(len(level_vector))])


def get_num_dof_of_full_grids(level_vectors, boundary) -> int:
    num_dof = 0
    for level_vector in level_vectors:
        num_dof += get_num_dof_of_full_grid(level_vector, boundary)
    return num_dof


def get_num_dof_of_subspace(level_vector, boundary) -> int:
    return np.prod([2**(level_vector[i]-1) if level_vector[i] > 0 else boundary[i] for i in range(len(level_vector))])


def get_num_dof_of_subspaces(level_vectors: set(tuple), boundary) -> int:
    num_dof = 0
    for level_vector in level_vectors:
        num_dof += get_num_dof_of_subspace(level_vector, boundary)
    return num_dof


# computes the active set by minimum level difference
# todo also implement more sensible schemes like the tilted plane by Christoph Kowitz
def compute_active_set(lmin, lmax):
    dim = len(lmin)
    firstLevelDifference = lmax[0] - lmin[0]
    uniformLevelDifference = [(lmax[i] - lmin[i]) ==
                              firstLevelDifference for i in range(dim)]
    diagonalIndex = 0
    s = set()
    if uniformLevelDifference and firstLevelDifference >= dim:
        levelSum = sum(lmin) + firstLevelDifference - diagonalIndex
        for offset in partition_integer_in_num_partitions_with_zeros(firstLevelDifference-diagonalIndex, dim):
            # ic(offset)
            grid = np.array(lmin) + np.array(offset)
            assert (np.sum(grid) == levelSum)
            if (np.array(grid) >= np.array(lmin)).all():
                s.add(tuple(grid))
    else:
        listOfRanges = [list(range(0, lmax[i]+1))
                        for i in range(len(lmax))]
        listOfAllGridsUpToLmax = list(it.product(*listOfRanges))
        levelSum = get_min_level_sum(lmin, lmax)+diagonalIndex
        for grid in listOfAllGridsUpToLmax:
            if (np.sum(grid) == levelSum and (np.array(grid) >= np.array(lmin)).all()):
                s.add(grid)

    return s


def get_graph_from_active_set(main_diagonal) -> nx.Graph:
    """ get a networkx graph representation of the active set
        where the graph's key labels are strings of level vectors
    """
    G = nx.Graph()

    # the maximum possible number of connections is given by the binomial coefficient
    # (select two possible dimensions, twice for each sign of +-1 to keep the same level sum)
    max_number_connections = int(2*binom(len(main_diagonal[0]), 2))
    main_diagonal_list_of_lists = [list(level) for level in main_diagonal]
    for level in main_diagonal:
        # for every level on the main diagonal / in the active set, add a node to the graph
        G.add_node(str(level), node_size=1)
        # and add an edge between all nodes where the level vectors differ by 2
        # only works for "regular" diagonal active sets
        offset = np.array([0]*len(level), dtype=int)
        offset[0] = 1
        offset[1] = -1
        for permutation in it.permutations(offset):
            level2 = level+permutation
            assert int(np.linalg.norm(level-level2, ord=1)) == 2
            if level2.tolist() in main_diagonal_list_of_lists:
                G.add_edge(str(level), str(level2), edge_value=1, capacity=1)
                # ic("added edge", level, level2)
        num_connections = len(G.edges(str(level)))
        assert (num_connections > 0)
        assert (num_connections <= max_number_connections)

    assert (nx.is_connected(G))
    return G


def compute_regular_combination_dictionary(lmin, levelDifference: int) -> dict:
    # implements the standard formula, cf. "Sparse Grids in a Nutshell"
    # (cf. https://link.springer.com/chapter/10.1007/978-3-642-31703-3_3)
    dim = len(lmin)
    combination_dictionary = {}
    for q in range(dim):
        coeff = (-1)**q * binom(dim-1, q)
        levelSum = sum(lmin) + levelDifference - q
        # ic(coeff, levelSum)
        for offset in partition_integer_in_num_partitions_with_zeros(levelDifference-q, dim):
            # ic(offset)
            grid = np.array(lmin) + np.array(offset)
            if (np.sum(grid) == levelSum and (np.array(grid) >= np.array(lmin)).all()):
                combination_dictionary[tuple(grid)] = coeff
    return combination_dictionary


def compute_adaptive_combination_dictionary(active_set) -> dict:
    # algorithm that can be derived by Hamming distance
    # (cf. Harding 2016, https://link.springer.com/content/pdf/10.1007%2F978-3-319-28262-6_4.pdf)
    combination_dictionary = {}
    for l in active_set:
        combination_dictionary[l] = 1

    dict_of_subspaces = {}
    for l in combination_dictionary:
        for subspace in get_downward_closed_set_from_level_vector(l):
            if subspace in dict_of_subspaces:
                dict_of_subspaces[subspace] += 1
            else:
                dict_of_subspaces[subspace] = 1

    # # remove subspaces which are too much
    while (set(dict_of_subspaces.values()) != set([1])):
        for subspace in dict_of_subspaces:
            currentCount = dict_of_subspaces[subspace]
            if currentCount != 1:
                diff = currentCount - 1

                if subspace in combination_dictionary:
                    combination_dictionary[subspace] -= diff
                    if combination_dictionary[subspace] == 0:
                        del combination_dictionary[subspace]
                else:
                    combination_dictionary[subspace] = -diff

                for l in get_downward_closed_set_from_level_vector(subspace):
                    dict_of_subspaces[l] -= diff
    return combination_dictionary


def compute_combination_dictionary(lmin, lmax) -> dict:
    dim = len(lmin)
    firstLevelDifference = lmax[0] - lmin[0]
    uniformLevelDifference = [(lmax[i] - lmin[i]) ==
                              firstLevelDifference for i in range(dim)]
    if uniformLevelDifference and firstLevelDifference >= dim:
        ic("binomial")
        return compute_regular_combination_dictionary(lmin, firstLevelDifference)

    else:
        ic("non-binomial")
        active_set = compute_active_set(lmin, lmax)
        ic("warning: no other active set implemented yet, will give same result as binomial but take longer!")
        return compute_adaptive_combination_dictionary(active_set)


def get_combination_dictionary_from_file(filename) -> dict:
    combination_dictionary = {}
    data = combischeme_output.read_data_from_json(filename)
    for grid in data:
        combination_dictionary[tuple(grid['level'])] = grid['coeff']
    return combination_dictionary


class CombinationScheme():
    def get_lmax(self) -> list[int]:
        return self._lmax

    def get_lmin(self) -> list[int]:
        return self._lmin

    def get_dimensionality(self) -> int:
        return len(self._lmin)

    def get_num_component_grids(self) -> int:
        return len(self._combination_dictionary)

    def get_boundary_points(self) -> list:
        assert self._boundary_points is not None
        return self._boundary_points

    def get_total_num_points_combi(self) -> int:
        total_num_points = 0
        for level in self.get_levels_of_nonzero_coefficient():
            total_num_points += get_num_dof_of_full_grid(
                level, self._boundary_points)
        return total_num_points

    def get_num_grids_per_level_sum(self) -> dict:
        num_grids_per_level_sum = {}
        for level in self.get_levels_of_nonzero_coefficient():
            levelSum = np.sum([l for l in level])
            if levelSum in num_grids_per_level_sum:
                num_grids_per_level_sum[levelSum] += 1
            else:
                num_grids_per_level_sum[levelSum] = 1
        return num_grids_per_level_sum

    def get_combination_dictionary(self) -> dict:
        return self._combination_dictionary

    def get_levels_of_nonzero_coefficient(self):
        return self._combination_dictionary.keys()

    def get_nonzero_coefficients(self):
        # the active front will have coefficients = 1
        # (but they may not be the only ones)
        return self._combination_dictionary.values()

    def get_coefficient(self, level: tuple) -> float:
        return self._combination_dictionary[level]

    def get_sparse_grid_spaces(self) -> set:
        spaces = get_downward_closed_set_from_level_vectors(
            self.get_levels_of_nonzero_coefficient())
        assert spaces is not None
        assert spaces != set()
        return spaces

    def get_necessary_sparse_grid_spaces(self) -> set:
        raise NotImplementedError(
            "get_necessary_sparse_grid_spaces was not correctly implemented for this class")

    def add_component(self, level: tuple, coefficient: float) -> None:
        self._combination_dictionary[level] = coefficient


class CombinationSchemeFromMaxLevel(CombinationScheme):
    def __init__(self, lmax, lmin=None, boundary_points=None):
        self._lmax = lmax
        if lmin is None:
            self._lmin = [1]*len(self.lmax)
        else:
            self._lmin = lmin
        if boundary_points is None:
            self._boundary_points = [2]*len(self._lmax)
        else:
            self._boundary_points = boundary_points
        assert (len(self._lmin) == len(self._lmax))
        for i in range(len(lmax)):
            assert (lmin[i] <= lmax[i])
        self._combination_dictionary = compute_combination_dictionary(
            lmin, lmax)
        assert (self._combination_dictionary is not None)
        self._graph = None

    def get_necessary_sparse_grid_spaces(self) -> set:
        lmax_reduced = [max(self._lmax[i]-1, self._lmin[i])
                        for i in range(len(self._lmax))]
        return get_downward_closed_set_from_level_vectors(compute_active_set(self._lmin, lmax_reduced))

    def assign_remaining_with_little_intersection_but_load_balanced(self, partitioned_schemes: list(CombinationScheme)) -> list(CombinationScheme):
        levels = self.get_levels_of_nonzero_coefficient()
        # filter out levels that are already assigned
        for scheme in partitioned_schemes:
            # currently I am assuming that only the active front is already assigned
            assert all([self.get_coefficient(
                l) == 1 for l in levels if l in scheme.get_levels_of_nonzero_coefficient()])
            assert all([self.get_coefficient(l) == scheme.get_coefficient(
                l) for l in levels if l in scheme.get_levels_of_nonzero_coefficient()])
            levels = [
                l for l in levels if l not in scheme.get_levels_of_nonzero_coefficient()]
        assert (len(levels) < len(
            self.get_levels_of_nonzero_coefficient()))

        # assign lower diagonals
        # count number of partial schemes that shadow the level
        shadow_dict = {}
        for scheme_index, scheme in enumerate(partitioned_schemes):
            part_levels = scheme.get_levels_of_nonzero_coefficient()
            part_subspaces = scheme.get_sparse_grid_spaces()
            assert len(part_levels) > 0
            assert len(part_subspaces) >= len(part_levels)
            for level in levels:
                if level not in shadow_dict:
                    shadow_dict[level] = []
                # keep track of which partitions level could be assigned to (at absolutely no added sparse grid cost)
                if level in part_subspaces:
                    shadow_dict[level].append(scheme_index)

        # assert that all levels can go somewhere for free
        assert all([(len(shadow_dict[level]) > 0) for level in shadow_dict])
        max_num_shadowing = max([len(shadow_dict[level])
                                for level in shadow_dict])
        single_shadowed = [
            level for level in shadow_dict if len(shadow_dict[level]) == 1]
        # the subspaces belonging to single-shadowed grids are interesting because they can be left out of the between-partition communication
        ic(len(single_shadowed), get_num_dof_of_subspaces(
            single_shadowed, self.get_boundary_points()))
        for level in single_shadowed:
            # assign to partition
            partitioned_schemes[shadow_dict[level][0]].add_component(
                level, self.get_coefficient(level))
            del shadow_dict[level]

        # from here on, it is heuristic
        # take level sum as proxy for size
        level_sizes = set([sum(l) for l in shadow_dict])
        levels_by_size = {size: [l for l in shadow_dict if sum(
            l) == size] for size in level_sizes}
        shadow_dict_copy = shadow_dict.copy()

        # start with the largest grids
        for level_size in sorted(level_sizes, reverse=True):
            for level in levels_by_size[level_size]:
                # go from those grids with few options to those with many options
                for num_shadowing in range(2, max_num_shadowing+1):
                    if len(shadow_dict_copy[level]) == num_shadowing:
                        # assign to partition with the lowest current number of DOF
                        candidates = shadow_dict[level]
                        candidates.sort(reverse=False, key=lambda i: get_num_dof_of_full_grids(
                            partitioned_schemes[i].get_levels_of_nonzero_coefficient(), self.get_boundary_points()))
                        assign_to = candidates[0]
                        partitioned_schemes[assign_to].add_component(
                            level, self.get_coefficient(level))
                        del shadow_dict[level]
        assert (len(shadow_dict) == 0)

        num_dof_per_partition = [get_num_dof_of_full_grids(scheme.get_levels_of_nonzero_coefficient(
        ), self.get_boundary_points()) for scheme in partitioned_schemes]
        num_grids_per_partition = [
            len(scheme.get_levels_of_nonzero_coefficient()) for scheme in partitioned_schemes]
        ic(min(num_dof_per_partition), max(num_dof_per_partition))
        ic(min(num_grids_per_partition), max(num_grids_per_partition))

        # TODO add CombinationSchemeFromCombinationSchemes and check if they are the same as self
        # for now, just make sure that all levels are assigned
        test_levels = self.get_levels_of_nonzero_coefficient()
        for scheme in partitioned_schemes:
            test_levels = [
                l for l in test_levels if l not in scheme.get_levels_of_nonzero_coefficient()]
        assert (len(test_levels) == 0)

        return partitioned_schemes

    def get_graph(self) -> nx.Graph:
        if self._graph is None:
            self._main_diagonal = [np.array(l, dtype=int)
                                   for l in compute_active_set(self._lmin, self._lmax)]
            self._graph = get_graph_from_active_set(
                self._main_diagonal)
        return self._graph, self._main_diagonal

    def split_scheme_metis(self, num_partitions: int = 2,
                           objtype='vol',
                           ctype='shem',
                           iptype='grow',
                           rtype='fm',
                           ncuts=100,
                           niter=1000,
                           ufactor=100,
                           minconn=True,
                           **metis_kwargs) -> list(CombinationScheme):
        # import metis only if required
        import metis

        G, main_diagonal = self.get_graph()

        max_number_connections = int(2*binom(len(main_diagonal[0]), 2))
        for level in main_diagonal:
            num_connections = len(G.edges(str(level)))
            # add costs for nodes at the exterior of the active set
            # this is where you would also enter a load model possibly
            # G.nodes[str(level)]['node_value'] = 1 + 1 * \
            #     (max_number_connections - num_connections)**1

        # cf. https://stackoverflow.com/a/50686793
        # which node attributes should metis use for partitioning
        G.graph['node_weight_attr'] = 'node_value'
        G.graph['node_size_attr'] = 'node_size'
        G.graph['edge_weight_attr'] = 'edge_value'

        # partitions from METIS
        (cut, parts) = metis.part_graph(G, num_partitions,
                                        contig=True,
                                        objtype=objtype,
                                        ctype=ctype,
                                        iptype=iptype,
                                        rtype=rtype,
                                        ncuts=ncuts,
                                        niter=niter,
                                        ufactor=ufactor,
                                        minconn=minconn,
                                        seed=1111,
                                        ** metis_kwargs
                                        )
        # make sure we get all the partitions we asked for
        assert min(parts) == 0
        for i in range(num_partitions):
            assert i in parts
        assert max(parts) == num_partitions-1
        for i, node in enumerate(G.nodes()):
            G.nodes[node]['parts'] = parts[i]

        # draw graph -- visible in jupyter notebook
        nx.draw(G, pos=nx.drawing.layout.spring_layout(G),
                with_labels=True, node_color=parts, vmin=-1)

        # transfer back to the ordering we had in the main diagonal
        parts = [G.nodes[str(level)]['parts'] for level in main_diagonal]

        partitioned_schemes = [CombinationSchemeFromCombinationDictionary({tuple(main_diagonal[i]): 1 for i in range(
            len(main_diagonal)) if parts[i] == p}, boundary_points=self.get_boundary_points()) for p in range(num_partitions)]

        partitioned_schemes = self.assign_remaining_with_little_intersection_but_load_balanced(
            partitioned_schemes)

        return partitioned_schemes

    def split_scheme_level_sum(self, num_partitions: int = 2) -> list(CombinationScheme):
        assert (num_partitions == 2)
        main_diagonal = [np.array(l, dtype=int)
                         for l in compute_active_set(self._lmin, self._lmax)]
        G = get_graph_from_active_set(main_diagonal)

        dim = len(main_diagonal[0])
        dims_system_1 = slice(0, dim, num_partitions)
        dims_system_2 = slice(1, dim, num_partitions)
        for level in main_diagonal:
            # assign by level sum
            lower_dims_diff_to_lmax = np.array(
                self._lmax[dims_system_1])-np.array(level)[dims_system_1]
            norm_lower_dims = - \
                np.linalg.norm(lower_dims_diff_to_lmax, ord=0.99)
            higher_dims_diff_to_lmax = np.array(
                self._lmax[dims_system_2])-np.array(level[dims_system_2])
            norm_higher_dims = - \
                np.linalg.norm(higher_dims_diff_to_lmax, ord=0.99)
            if norm_lower_dims > norm_higher_dims:
                G.nodes[str(level)]['parts'] = 0
            else:
                G.nodes[str(level)]['parts'] = 1
        # draw graph -- visible in jupyter notebook
        colors = [G.nodes[node]['parts'] for node in list(G.nodes())]
        nx.draw(G, pos=nx.drawing.layout.spring_layout(G),
                with_labels=True, node_color=colors, vmin=-1)

        parts = [G.nodes[str(level)]['parts'] for level in main_diagonal]
        assert min(parts) == 0
        assert max(parts) == num_partitions-1
        partitioned_schemes = [CombinationSchemeFromCombinationDictionary({tuple(main_diagonal[i]): 1 for i in range(
            len(main_diagonal)) if parts[i] == p}, boundary_points=self.get_boundary_points()) for p in range(num_partitions)]

        partitioned_schemes = self.assign_remaining_with_little_intersection_but_load_balanced(
            partitioned_schemes)

        return partitioned_schemes


class CombinationSchemeFromCombinationDictionary(CombinationScheme):
    def __init__(self, dictionary, boundary_points=None):
        self._combination_dictionary = dictionary
        assert (self._combination_dictionary is not None)
        self._lmax = np.max(
            list(self.get_levels_of_nonzero_coefficient()), axis=0)
        self._lmin = np.min(
            list(self.get_levels_of_nonzero_coefficient()), axis=0)
        ic(self._lmax, self._lmin)
        assert (len(self._lmin) == len(self._lmax))
        for i in range(len(self._lmax)):
            assert (self._lmin[i] <= self._lmax[i])
        if boundary_points is None:
            self._boundary_points = [2]*len(self._lmax)
        else:
            self._boundary_points = boundary_points


class CombinationSchemeFromFile(CombinationSchemeFromCombinationDictionary):
    def __init__(self, filename: str, boundary_points=None):
        dictionary = get_combination_dictionary_from_file(
            filename)
        # if boundary_points is a scalar, assume it is the same for all dimensions
        if isinstance(boundary_points, int):
            boundary_points = [boundary_points] * \
                len(next(iter(dictionary.keys())))
        super().__init__(dictionary, boundary_points)


def write_scheme_to_json(scheme: CombinationScheme, file_name: str = None):
    if file_name is None:
        # assuming double data type = 8 bytes = 64 bit
        mem = (scheme.get_total_num_points_combi()*8)
        # ic(readable_bytes(mem))
        file_name = 'scheme_' + \
            combischeme_output.readable_bytes(mem) + '.json'
    combischeme_output.write_scheme_dictionary_to_json(
        scheme.get_combination_dictionary(), file_name)


def get_conjoint_subspaces(schemes: list(CombinationScheme)) -> set:
    p_subspaces = [scheme.get_sparse_grid_spaces() for scheme in schemes]
    all_conjoint = set(reduce(and_, p_subspaces))
    return all_conjoint


def get_union_subspaces(schemes: list(CombinationScheme)) -> set:
    p_subspaces = [scheme.get_sparse_grid_spaces() for scheme in schemes]
    all_union = set(reduce(or_, p_subspaces))
    return all_union


def split_scheme_by_level_sum(scheme: CombinationScheme) -> tuple(CombinationScheme, CombinationScheme):
    """splits a combination scheme into two: input is assumed to be valid regular combination scheme"""
    gridsForSystems = list(scheme.get_levels_of_nonzero_coefficient())
    dim = scheme.get_dimensionality()
    lmax = scheme.get_lmax()
    gridsForSystem1 = []
    gridsForSystem2 = []
    ic(len(gridsForSystem1), len(gridsForSystem2), len(gridsForSystems))

    gridsToIterate = gridsForSystems.copy()
    rr = 0
    # # for different splits -- the sum of weights is assumed to be 1.,
    # # if system 1's weight is 0.5, we have an even split.
    weight_system_1 = 0.5
    dims_system_1 = slice(0, dim, 2)
    dims_system_2 = slice(1, dim, 2)
    for level in gridsToIterate:
        lower_dims_diff_to_lmax = np.array(
            lmax[dims_system_1])-np.array(level)[dims_system_1]
        norm_lower_dims = - \
            np.linalg.norm(lower_dims_diff_to_lmax, ord=0.99)
        higher_dims_diff_to_lmax = np.array(
            lmax[dims_system_2])-np.array(level[dims_system_2])
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

    assert (len(gridsForSystems) == 0)
    ic(len(gridsForSystem1), len(gridsForSystem2), len(gridsForSystems))

    # build the new combination dictionaries
    dictionary1 = {}
    for level in gridsForSystem1:
        dictionary1[level] = scheme.get_coefficient(level)
    dictionary2 = {}
    for level in gridsForSystem2:
        dictionary2[level] = scheme.get_coefficient(level)
    return (CombinationSchemeFromCombinationDictionary(dictionary1, boundary_points=scheme.get_boundary_points()), CombinationSchemeFromCombinationDictionary(dictionary2, boundary_points=scheme.get_boundary_points()))


def split_scheme_by_single_dimension(scheme: CombinationScheme, from_this_level_and_above, dimension_to_split_at=-1) -> tuple(CombinationScheme, CombinationScheme):
    """splits a combination scheme into two: input is assumed to be valid regular combination scheme"""
    gridsForSystems = list(scheme.get_levels_of_nonzero_coefficient())
    gridsForSystem1 = []
    gridsForSystem2 = []
    ic(len(gridsForSystem1), len(gridsForSystem2), len(gridsForSystems))

    gridsToIterate = gridsForSystems.copy()
    for level in gridsToIterate:
        if level[dimension_to_split_at] >= from_this_level_and_above:
            gridsForSystems.remove(level)
            gridsForSystem1.append(level)
            continue
        else:
            gridsForSystems.remove(level)
            gridsForSystem2.append(level)

    assert (len(gridsForSystems) == 0)
    ic(len(gridsForSystem1), len(gridsForSystem2), len(gridsForSystems))

    # build the new combination dictionaries
    dictionary1 = {}
    for level in gridsForSystem1:
        dictionary1[level] = scheme.get_coefficient(level)
    dictionary2 = {}
    for level in gridsForSystem2:
        dictionary2[level] = scheme.get_coefficient(level)
    return (CombinationSchemeFromCombinationDictionary(dictionary1, boundary_points=scheme.get_boundary_points()), CombinationSchemeFromCombinationDictionary(dictionary2, boundary_points=scheme.get_boundary_points()))


def assign_combischeme_to_groups(scheme: CombinationScheme, num_process_groups: int, num_dof_buffer_for_first_group=0) -> list(dict):
    dim = scheme.get_dimensionality()
    total_num_points_combi = scheme.get_total_num_points_combi()

    levels = list(scheme.get_levels_of_nonzero_coefficient())
    levels.sort(key=lambda x: get_num_dof_of_full_grid(
        x, [2]*dim), reverse=True)
    ic(levels[:5])

    assignment = []
    for i in range(num_process_groups):
        assignment.append({})
    assigned_FG_size = [0]*num_process_groups
    assigned_FG_size[0] = num_dof_buffer_for_first_group

    nextIndex = num_process_groups-1
    for level in levels:
        assignment[nextIndex][level] = scheme.get_coefficient(level)
        assigned_FG_size[nextIndex] += get_num_dof_of_full_grid(
            level, scheme.get_boundary_points())
        # this is where load balancing happens!
        nextIndex = np.argmin(assigned_FG_size)

    assigned_FG_size[0] -= num_dof_buffer_for_first_group
    ic(assigned_FG_size[0], max(assigned_FG_size), min(assigned_FG_size))
    ic(max([len(a) for a in assignment]), min([len(a) for a in assignment]))
    assert (sum(assigned_FG_size) == total_num_points_combi)
    assert (sum([len(a) for a in assignment]) == len(levels))
    return assignment, assigned_FG_size
