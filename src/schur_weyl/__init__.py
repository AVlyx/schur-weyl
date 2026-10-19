r"""Partitions, Young tableaux, symmetric-group characters and Schur–Weyl duality,
built from scratch on `numpy` and the standard library — no
representation-theory or combinatorics dependencies.

Schur–Weyl duality splits $(\mathbb{C}^d)^{\otimes n}$ under the commuting
actions of $S_n$ (permuting the tensor factors) and $GL_d$ (acting diagonally)
into blocks indexed by partitions $\lambda \vdash n$ with $\ell(\lambda) \le d$:

$$
(\mathbb{C}^d)^{\otimes n} \;\cong\; \bigoplus_{\lambda} S^\lambda \otimes V_\lambda^d .
$$

This package gives you both sides of that sum — the combinatorics of $\lambda$,
and the linear algebra that projects onto each block.

# Install

```bash
pip install schur-weyl
```

Requires Python ≥ 3.9. The only runtime dependency is `numpy`. For local
development, install it editable with the `dev` extra, which adds `pytest` and
`ruff`:

```bash
pip install -e ".[dev]"
```

# Quick start

Partitions are plain tuples, tableaux are tuples of tuples, and permutations are
0-indexed one-line tuples. Nothing is a custom class, so everything composes
with ordinary Python.

```python
from schur_weyl import partitions, dim_specht, dim_weyl, isotypic_proj
from schur_weyl.character import character
from schur_weyl.sw_measure import schur_weyl_measure

[*partitions(4)]
# [(4,), (3, 1), (2, 2), (2, 1, 1), (1, 1, 1, 1)]

dim_specht((2, 1))            # 2   -- the standard rep of S_3
dim_weyl((2, 1), d=3)         # 8

character((2, 1), (3,))       # -1  -- Murnaghan-Nakayama

schur_weyl_measure([0.7, 0.3], k=4)
# {(4,): 0.4141, (3, 1): 0.4977, (2, 2): 0.0882, ...}

P = isotypic_proj((2, 1), d=2)   # np.ndarray, P @ P == P
```

Every function carries a docstring with worked examples, run as doctests. The
test suite checks closed-form identities — $\sum_\lambda (f^\lambda)^2 = n!$,
character orthogonality, $P_\lambda^2 = P_\lambda$ — rather than fixed outputs.

`schur-weyl` is developed at
[AVlyx/schur-weyl](https://github.com/AVlyx/schur-weyl), released on
[PyPI](https://pypi.org/project/schur-weyl/), and licensed under the MIT
License.
"""

from .young_diagrams import partitions, partition_conjugate, hook_lenght, cells, majorizes, addable_corners, removable_corners, is_in_partition
from .dimensions import dim_specht, dim_weyl
from .tableaux import standard_young_tableaux, semi_standard_young_tableau_by_content, semi_standard_young_tableau, kostka, kostka_matrix, inverse_kostka_matrix, reverse_reading_word
from .symmetric_functions import schur_polynomial, power_sum, monomial_symmetric, monomial_in_schur_basis
from .isotypic import isotypic_proj, apply_isotypic_proj
from .kronecker import kronecker_coefficient
from .schur_transform import schur_transform, schur_transform_matrix
