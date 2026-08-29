"""The roadmap gates: the identities that must hold exactly, for the phases implemented so far."""

from math import factorial

import pytest

from schur_weyl.dimensions import dim_specht, dim_weyl
from schur_weyl.tableaux import standard_young_tableaux
from schur_weyl.young_diagrams import partition_conjugate, partitions

PARTITION_COUNTS = [1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42]  # n = 0..10


# Phase 0 -- partitions and Young diagrams


@pytest.mark.parametrize("n, count", list(enumerate(PARTITION_COUNTS)))
def test_partition_count(n, count):
    assert len(list(partitions(n))) == count


@pytest.mark.parametrize("n", range(11))
def test_conjugate_is_an_involution(n):
    for lam in partitions(n):
        assert partition_conjugate(partition_conjugate(lam)) == lam


# Phase 1 -- dimension formulas


@pytest.mark.parametrize("n", range(10))
def test_sum_of_squared_specht_dims_is_n_factorial(n):
    assert sum(dim_specht(lam) ** 2 for lam in partitions(n)) == factorial(n)


@pytest.mark.parametrize("d", [2, 3, 4, 5])
@pytest.mark.parametrize("n", range(10))
def test_schur_weyl_dimension_count(n, d):
    assert sum(dim_specht(lam) * dim_weyl(lam, d) for lam in partitions(n, max_height=d)) == d**n


# Phase 2 -- tableaux


@pytest.mark.parametrize("n", range(8))
def test_syt_count_matches_specht_dim(n):
    for lam in partitions(n):
        assert len(list(standard_young_tableaux(lam))) == dim_specht(lam)
