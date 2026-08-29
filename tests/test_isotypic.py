from functools import reduce
from itertools import combinations
from math import prod

import numpy as np
import pytest

from schur_weyl.dimensions import dim_specht, dim_weyl
from schur_weyl.isotypic import _perm_index_map, isotypic_proj
from schur_weyl.symmetric_group import all_permutations, permutation_compose, permutation_cycle_type
from schur_weyl.young_diagrams import partitions


def perm_action(psi: np.ndarray, sigma: tuple[int, ...], d: int) -> np.ndarray:
    """P(sigma) applied to a flattened state of (C^d)^{otimes n}."""
    return psi[_perm_index_map(sigma, d)]


# Phase 6 -- tensor space and the permutation action


@pytest.mark.parametrize("d", [2, 3])
def test_permutation_action_is_a_homomorphism(d):
    n = 4
    rng = np.random.default_rng(0)
    psi = rng.normal(size=d**n) + 1j * rng.normal(size=d**n)
    for sigma in all_permutations(n):
        for tau in all_permutations(n):
            composed = perm_action(psi, permutation_compose(sigma, tau), d)
            assert np.array_equal(perm_action(perm_action(psi, tau, d), sigma, d), composed)


def test_permuted_tensor_power_trace():
    n, d = 4, 3
    rng = np.random.default_rng(0)
    a = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    rho = a @ a.conj().T
    rho /= np.trace(rho).real
    rho_n = reduce(np.kron, [rho] * n)
    for sigma in all_permutations(n):
        trace = rho_n[_perm_index_map(sigma, d), np.arange(d**n)].sum()
        cycles = prod(np.trace(np.linalg.matrix_power(rho, c)) for c in permutation_cycle_type(sigma))
        assert np.isclose(trace, cycles)


# Phase 7 -- isotypic projectors


@pytest.mark.parametrize("n, d", [(4, 2), (4, 3), (5, 3)])
def test_isotypic_projectors(n, d):
    projectors = {lam: isotypic_proj(lam, d) for lam in partitions(n)}
    total = np.zeros((d**n, d**n))
    for lam, p in projectors.items():
        assert np.allclose(p @ p, p)
        assert np.isclose(np.trace(p), dim_specht(lam) * dim_weyl(lam, d))
        total += p
    assert np.allclose(total, np.eye(d**n))
    for lam, mu in combinations(projectors, 2):
        assert np.allclose(projectors[lam] @ projectors[mu], 0)
