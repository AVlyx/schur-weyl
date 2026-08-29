from math import factorial

import pytest

from schur_weyl.symmetric_group import all_permutations_by_cycle_type, permutation_conjugacy_class_size
from schur_weyl.young_diagrams import partitions


@pytest.mark.parametrize("n", range(1, 9))
def test_conjugacy_class_sizes_sum_to_n_factorial(n):
    assert sum(permutation_conjugacy_class_size(mu) for mu in partitions(n)) == factorial(n)


@pytest.mark.parametrize("n", range(1, 9))
def test_conjugacy_class_sizes_match_brute_force(n):
    classes = all_permutations_by_cycle_type(n)
    for mu in partitions(n):
        assert len(classes[mu]) == permutation_conjugacy_class_size(mu)
