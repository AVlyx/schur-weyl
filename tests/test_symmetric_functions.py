import numpy as np
import pytest

from schur_weyl.character import character
from schur_weyl.dimensions import dim_weyl
from schur_weyl.symmetric_functions import power_sum, schur_polynomial
from schur_weyl.young_diagrams import partitions


@pytest.mark.parametrize("n", range(1, 7))
def test_power_sum_expands_over_characters(n):
    xs = list(np.random.default_rng(0).random(4))
    for mu in partitions(n):
        expansion = sum(character(lam, mu) * schur_polynomial(lam, xs) for lam in partitions(n))
        assert np.isclose(power_sum(mu, xs), expansion)


@pytest.mark.parametrize("d", [2, 3, 4])
@pytest.mark.parametrize("n", range(1, 7))
def test_schur_at_all_ones_is_weyl_dim(n, d):
    for lam in partitions(n):
        assert np.isclose(schur_polynomial(lam, [1.0] * d), dim_weyl(lam, d))
