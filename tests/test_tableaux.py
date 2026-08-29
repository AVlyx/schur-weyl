from math import factorial, prod

import pytest

from schur_weyl.dimensions import dim_specht, dim_weyl
from schur_weyl.tableaux import kostka, semi_standard_young_tableau, standard_young_tableaux
from schur_weyl.young_diagrams import majorizes, partitions


@pytest.mark.parametrize("n", range(8))
def test_syt_count_matches_specht_dim(n):
    for lam in partitions(n):
        assert len(list(standard_young_tableaux(lam))) == dim_specht(lam)


@pytest.mark.parametrize("d", [2, 3])
@pytest.mark.parametrize("n", range(7))
def test_ssyt_count_matches_weyl_dim(n, d):
    for lam in partitions(n):
        assert len(list(semi_standard_young_tableau(lam, d))) == dim_weyl(lam, d)


@pytest.mark.parametrize("n", range(1, 7))
def test_kostka_vs_dominance(n):
    for lam in partitions(n):
        assert kostka(lam, (1,) * n) == dim_specht(lam)
        assert kostka(lam, lam) == 1
        for mu in partitions(n):
            assert (kostka(lam, mu) != 0) == majorizes(lam, mu)


@pytest.mark.parametrize("n", range(1, 7))
def test_kostka_weighted_sum_is_multinomial(n):
    for mu in partitions(n):
        multinomial = factorial(n) // prod(factorial(m) for m in mu)
        assert sum(kostka(lam, mu) * dim_specht(lam) for lam in partitions(n)) == multinomial
