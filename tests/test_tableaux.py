from itertools import permutations, product
from math import factorial, prod

import numpy as np
import pytest

from schur_weyl.dimensions import dim_specht, dim_weyl
from schur_weyl.tableaux import _unitriangular_inverse, inverse_kostka_matrix, kostka, kostka_matrix, semi_standard_young_tableau, semi_standard_young_tableau_by_content, standard_young_tableaux
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


def _brute_force_kostka(lam, mu):
    return len(list(semi_standard_young_tableau_by_content(lam, mu)))


@pytest.mark.parametrize("n", range(8))
def test_kostka_matches_brute_force(n):
    for lam in partitions(n):
        for mu in partitions(n):
            k = kostka(lam, mu)
            assert type(k) is int
            assert k == _brute_force_kostka(lam, mu)


@pytest.mark.parametrize("n", range(1, 6))
def test_kostka_compositions_match_brute_force(n):
    # every weak composition of n into 4 parts: unsorted contents with zero entries
    weights = [mu for mu in product(range(n + 1), repeat=4) if sum(mu) == n]
    for lam in partitions(n):
        for mu in weights:
            k = kostka(lam, mu)
            assert type(k) is int
            assert k == _brute_force_kostka(lam, mu)


@pytest.mark.parametrize("n", range(1, 7))
def test_kostka_symmetric_in_content(n):
    for lam in partitions(n):
        for mu in partitions(n):
            k = kostka(lam, mu)
            assert all(kostka(lam, p) == k for p in set(permutations(mu)))


@pytest.mark.parametrize("n", range(1, 21))
def test_kostka_standard_content_is_specht_dim(n):
    for lam in partitions(n, max_height=4):
        assert kostka(lam, (1,) * n) == dim_specht(lam)


@pytest.mark.parametrize("d", [2, 3, 4])
@pytest.mark.parametrize("n", range(1, 9))
def test_kostka_sum_over_weights_is_weyl_dim(n, d):
    # s_lam(1^d) = sum over weak compositions mu of n into d parts of K_{lam, mu}
    weights = [mu for mu in product(range(n + 1), repeat=d) if sum(mu) == n]
    for lam in partitions(n, max_height=d):
        assert sum(kostka(lam, mu) for mu in weights) == dim_weyl(lam, d)


def test_kostka_edge_cases():
    assert kostka((), ()) == 1
    assert kostka((), (0, 0)) == 1
    assert kostka((2, 1), (2,)) == 0  # sizes differ
    assert kostka((1, 1, 1), (2, 1)) == 0  # too many rows
    assert kostka((2, 1, 0), (1, 1, 1)) == 2
    assert type(kostka((3,), (1, 2))) is int


def test_kostka_large():
    lam = (8, 6, 4, 2)
    assert kostka(lam, (1,) * 20) == dim_specht(lam)
    assert kostka((12,), (3, 4, 5)) == 1


@pytest.mark.parametrize("n", range(10))
def test_kostka_matrix_matches_kostka(n):
    ps, table = kostka_matrix(n)
    assert ps == tuple(partitions(n))
    for a, lam in enumerate(ps):
        for b, mu in enumerate(ps):
            assert type(table[a][b]) is int
            assert table[a][b] == kostka(lam, mu)


@pytest.mark.parametrize("n", range(1, 13))
def test_kostka_matrix_unitriangular(n):
    # partitions(n) is a linear extension of dominance, so the matrix is upper unitriangular
    ps, table = kostka_matrix(n)
    for a in range(len(ps)):
        assert table[a][a] == 1
        assert all(table[a][b] == 0 for b in range(a))
        assert table[a][-1] == dim_specht(ps[a])


@pytest.mark.parametrize("n", range(13))
def test_inverse_kostka_matrix_is_inverse(n):
    ps, table = kostka_matrix(n)
    ps_inv, inverse = inverse_kostka_matrix(n)
    assert ps_inv == ps
    p = len(ps)
    for a in range(p):
        assert all(type(entry) is int for entry in inverse[a])
        for c in range(p):
            assert sum(inverse[a][b] * table[b][c] for b in range(p)) == (a == c)


@pytest.mark.parametrize("n", range(1, 11))
def test_unitriangular_inverse_python_ints_match_int64(n):
    _, table = kostka_matrix(n)
    assert _unitriangular_inverse(table, object).tolist() == _unitriangular_inverse(table, np.int64).tolist()


def test_unitriangular_inverse_detects_int64_overflow():
    big = 2**40
    table = ((1, big, 0), (0, 1, big), (0, 0, 1))  # inverse has the entry big**2 = 2**80
    assert _unitriangular_inverse(table, np.int64) is None
    assert _unitriangular_inverse(table, object).tolist() == [[1, -big, big**2], [0, 1, -big], [0, 0, 1]]
    too_big = ((1, 2**70), (0, 1))  # does not even fit in int64
    assert _unitriangular_inverse(too_big, np.int64) is None
