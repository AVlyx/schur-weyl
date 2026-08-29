"""Phase 4 gates -- characters via Murnaghan-Nakayama."""

from math import factorial

import pytest

from schur_weyl.character import character
from schur_weyl.dimensions import dim_specht
from schur_weyl.symmetric_group import permutation_conjugacy_class_size
from schur_weyl.young_diagrams import partition_conjugate, partitions


def sgn(mu: tuple[int, ...]) -> int:
    """The sign of a permutation of cycle type mu."""
    return (-1) ** (sum(mu) - len(mu))


@pytest.mark.parametrize("n", range(1, 9))
def test_character_at_identity_is_specht_dim(n):
    for lam in partitions(n):
        assert character(lam, (1,) * n) == dim_specht(lam)


@pytest.mark.parametrize("n", range(1, 9))
def test_trivial_and_sign_characters(n):
    for mu in partitions(n):
        assert character((n,), mu) == 1
        assert character((1,) * n, mu) == sgn(mu)


@pytest.mark.parametrize("n", range(1, 9))
def test_conjugate_character_is_sign_twisted(n):
    for lam in partitions(n):
        for mu in partitions(n):
            assert character(partition_conjugate(lam), mu) == sgn(mu) * character(lam, mu)


@pytest.mark.parametrize("n", range(1, 9))
def test_character_orthogonality(n):
    ps = list(partitions(n))
    sizes = [permutation_conjugacy_class_size(mu) for mu in ps]
    chi = [[character(lam, mu) for mu in ps] for lam in ps]
    for a in range(len(ps)):
        for b in range(len(ps)):
            # rows: sum_mu chi^a(mu) chi^b(mu) / z_mu == delta, cleared of denominators
            expected = factorial(n) if a == b else 0
            assert sum(sizes[c] * chi[a][c] * chi[b][c] for c in range(len(ps))) == expected
            # columns: sum_lam chi^lam(a) chi^lam(b) == delta * z_a
            expected = factorial(n) // sizes[a] if a == b else 0
            assert sum(chi[c][a] * chi[c][b] for c in range(len(ps))) == expected
