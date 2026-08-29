import numpy as np
import pytest

from schur_weyl.character import character
from schur_weyl.dimensions import dim_specht
from schur_weyl.symmetric_group import all_permutations, permutation_cycle_type
from schur_weyl.young_diagrams import partitions
from schur_weyl.young_orthonormal import (
    jucys_murphy,
    jucys_murphy_contents,
    young_orthogonal,
    young_orthogonal_generator,
)


@pytest.mark.parametrize("n", [4, 5])
def test_young_orthogonal_generator_relations(n):
    for lam in partitions(n):
        eye = np.eye(dim_specht(lam))
        s = [young_orthogonal_generator(lam, k) for k in range(1, n)]
        for k in range(n - 1):
            assert np.allclose(s[k] @ s[k].T, eye)
            assert np.allclose(s[k] @ s[k], eye)
        for k in range(n - 2):
            assert np.allclose(s[k] @ s[k + 1] @ s[k], s[k + 1] @ s[k] @ s[k + 1])
        for k in range(n - 1):
            for j in range(1, k - 1):
                assert np.allclose(s[j] @ s[k], s[k] @ s[j])


@pytest.mark.parametrize("n", [4, 5])
def test_young_orthogonal_traces_are_characters(n):
    for lam in partitions(n):
        for sigma in all_permutations(n):
            trace = np.trace(young_orthogonal(lam, sigma))
            assert np.isclose(trace, character(lam, permutation_cycle_type(sigma)))


@pytest.mark.parametrize("n", [4, 5])
def test_jucys_murphy_is_diagonal_of_contents(n):
    for lam in partitions(n):
        for k in range(1, n + 1):
            assert np.allclose(jucys_murphy(lam, k), np.diag(jucys_murphy_contents(lam, k)))
