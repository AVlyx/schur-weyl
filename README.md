# schur-weyl

Tools from Schur–Weyl duality, built from scratch on top of `numpy` and the Python standard library — no representation-theory or combinatorics dependencies.

Schur–Weyl duality decomposes $(\mathbb{C}^d)^{\otimes n}$ under the commuting actions of $S_n$ (permuting tensor factors) and $GL_d$ (acting diagonally) into isotypic blocks indexed by partitions $\lambda \vdash n$ with $\ell(\lambda) \le d$:

$$(\mathbb{C}^d)^{\otimes n} \cong \bigoplus_{\lambda} S^\lambda \otimes V_\lambda^d$$

This package implements the combinatorics and linear algebra needed to work with that decomposition directly: partitions and Young tableaux, dimension formulas, symmetric-group characters via Murnaghan–Nakayama, Schur polynomials, isotypic projectors on tensor space, the Schur–Weyl measure, and Young's orthogonal form for the irreducible representations of $S_n$.

## Installation

Requires Python ≥ 3.9.

```bash
pip install schur-weyl
```

The package itself only depends on `numpy`. For local development (cloning the repo, running the test suite), install it editable with the `dev` extra instead, which pulls in `pytest` and `ruff`:

```bash
pip install -e ".[dev]"
```

## Modules

| Module | Contents |
|---|---|
| [`young_diagrams`](src/schur_weyl/young_diagrams.py) | Partitions, conjugation, hook lengths, content, dominance order, addable/removable corners |
| [`dimensions`](src/schur_weyl/dimensions.py) | $f^\lambda$ (Specht module dimension, hook-length formula) and $\dim V_\lambda^d$ (Weyl module dimension, content formula) |
| [`tableaux`](src/schur_weyl/tableaux.py) | Standard and semistandard Young tableaux, Kostka numbers, reading words |
| [`symmetric_group`](src/schur_weyl/symmetric_group.py) | Permutations as 0-indexed one-line tuples: composition, inverse, cycle type, conjugacy class sizes |
| [`character`](src/schur_weyl/character.py) | Irreducible characters $\chi^\lambda(\mu)$ of $S_n$ via the Murnaghan–Nakayama rule; full character tables |
| [`symmetric_functions`](src/schur_weyl/symmetric_functions.py) | Schur polynomials (Jacobi–Trudi) and power-sum symmetric functions |
| [`isotypic`](src/schur_weyl/isotypic.py) | Isotypic projectors $P_\lambda$ acting on $(\mathbb{C}^d)^{\otimes n}$, built from class sums of the permutation action |
| [`sw_measure`](src/schur_weyl/sw_measure.py) | The Schur–Weyl measure $\Pr[\lambda] = f^\lambda\, s_\lambda(\mathrm{spec}\,\rho)$ for i.i.d. copies of a state $\rho$ |
| [`young_orthonormal`](src/schur_weyl/young_orthonormal.py) | Young's orthogonal form: explicit irrep matrices $\rho^\lambda(\sigma)$ indexed by standard tableaux, and the Jucys–Murphy elements |

## Usage

```python
from schur_weyl import (
    partitions, dim_specht, dim_weyl, kostka,
    schur_polynomial, isotypic_proj,
)
from schur_weyl.character import character, character_table
from schur_weyl.sw_measure import schur_weyl_measure

# dimension formulas
dim_specht((2, 1))          # -> 2, the standard rep of S_3
dim_weyl((2, 1), d=3)       # -> 8

# characters via Murnaghan-Nakayama
character((2, 1), (3,))     # -> -1

# the Schur-Weyl measure for n i.i.d. copies of a state with given spectrum
schur_weyl_measure(spectrum=[0.7, 0.3], k=4)
# -> {(4,): 0.4141, (3, 1): 0.4977, (2, 2): 0.0882, (2, 1, 1): 0.0, (1, 1, 1, 1): 0.0}

# the isotypic projector on (C^d)^{\otimes n}
P = isotypic_proj((2, 1), d=2)   # np.ndarray, P @ P == P
```

## Testing

```bash
pytest
```

Tests live in [`tests/`](tests) and mirror the module layout; most assert the closed-form identities from the roadmap (e.g. $\sum_\lambda (f^\lambda)^2 = n!$, character-table orthogonality, $P_\lambda^2 = P_\lambda$, $\operatorname{tr}\rho^\lambda(X_k)$ diagonal with content eigenvalues) rather than fixed expected outputs.

## Status

Implemented through Young's orthogonal form and Jucys–Murphy elements (roadmap phases 0–5, 7–9). Not yet implemented: the generic tensor-space permutation action as a standalone module (phase 6, currently inlined in `isotypic`), and the explicit Schur transform (phase 10). See [`TODO.md`](TODO.md) for near-term items.

## License

MIT — see [LICENSE](LICENSE).
