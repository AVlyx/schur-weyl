from functools import reduce

import numpy as np

from schur_weyl.isotypic import isotypic_proj
from schur_weyl.sw_measure import schur_weyl_measure


def test_schur_weyl_measure_matches_projector_traces():
    n, d = 4, 3
    spectrum = list(np.random.default_rng(0).dirichlet([1] * d))
    rho_n = reduce(np.kron, [np.diag(spectrum)] * n)
    probabilities = schur_weyl_measure(spectrum, n)
    assert np.isclose(sum(probabilities.values()), 1.0)
    for lam, p in probabilities.items():
        assert np.isclose(np.trace(isotypic_proj(lam, d) @ rho_n), p)
