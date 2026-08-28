"""Young's orthogonal form: explicit matrices for the irreps of S_n.

Conventions (consistent with the rest of the project):
  - cells are 0-indexed, content(i, j) = j - i
  - tableau ENTRIES are 1-indexed, i.e. 1..n
  - permutations are 0-indexed one-line tuples, compose(s, t)(i) = s(t(i))
  - the adjacent transposition s_k swaps the ENTRIES k and k+1, for
    1 <= k <= n-1.  It equals the 0-indexed permutation swapping k-1 and k.
"""

# Claude

from functools import lru_cache

import numpy as np
from tableaux import standard_young_tableaux
from dimensions import dim_specht

# ------------------------------------------------------------------ tableaux


@lru_cache(maxsize=None)
def syt_basis(lam: tuple[int, ...]) -> tuple[tuple[tuple[int, ...], ...], ...]:
    lam = tuple(x for x in lam if x > 0)
    return tuple(standard_young_tableaux(lam))


def _positions(T):
    """{entry: (i, j)} for a tableau."""
    return {v: (i, j) for i, row in enumerate(T) for j, v in enumerate(row)}


def _swap_entries(T, k):
    """The tableau with entries k and k+1 exchanged (may not be standard)."""
    rows = [list(r) for r in T]
    for r in rows:
        for c, v in enumerate(r):
            if v == k:
                r[c] = k + 1
            elif v == k + 1:
                r[c] = k
    return tuple(tuple(r) for r in rows)


# ---------------------------------------------------------- adjacent generator


def axial_distance(T, k):
    """cont(box holding k+1) - cont(box holding k), with content = col - row."""
    pos = _positions(T)
    i1, j1 = pos[k]
    i2, j2 = pos[k + 1]
    return (j2 - i2) - (j1 - i1)


@lru_cache(maxsize=None)
def young_orthogonal_generator(lam, k):
    """The matrix of rho^lam(s_k), where s_k swaps the entries k and k+1.

    On the span of {T, s_k T}, with r = 1 / axial_distance(T, k):

        [  r              sqrt(1 - r^2) ]
        [  sqrt(1 - r^2)  -r            ]

    If k and k+1 share a row the distance is +1 (so r = 1, eigenvalue +1);
    if they share a column it is -1 (eigenvalue -1).  In both cases s_k T is
    not standard and the off-diagonal term vanishes on its own.
    """
    lam = tuple(x for x in lam if x > 0)
    basis = syt_basis(lam)
    index = {T: a for a, T in enumerate(basis)}
    f = len(basis)
    M = np.zeros((f, f))

    for b, T in enumerate(basis):
        delta = axial_distance(T, k)
        r = 1.0 / delta
        M[b, b] = r
        if abs(delta) == 1:
            # same row (+1) or same column (-1): no mixing partner
            assert _swap_entries(T, k) not in index
            continue
        a = index[_swap_entries(T, k)]
        M[a, b] = np.sqrt(1.0 - r * r)
    return M


# ------------------------------------------------------------ general element


def _adjacent_factorization(sigma):
    """[i_1, ..., i_m] (0-indexed) with sigma == s_{i_1} s_{i_2} ... s_{i_m}."""
    arr = list(sigma)
    n = len(arr)
    right = []
    for i in range(n):
        for j in range(n - 1 - i):
            if arr[j] > arr[j + 1]:
                arr[j], arr[j + 1] = arr[j + 1], arr[j]
                right.append(j)
    return list(reversed(right))


def young_orthogonal(lam: tuple[int, ...], sigma: tuple[int, ...]):
    """rho^lam(sigma) for a 0-indexed one-line permutation tuple sigma.

    Factors sigma into adjacent transpositions and multiplies the generators.
    """
    lam = tuple(x for x in lam if x > 0)
    f = dim_specht(lam)
    M = np.eye(f)
    for i in _adjacent_factorization(sigma):
        M = M @ young_orthogonal_generator(lam, i + 1)  # 0-indexed i -> entry i+1
    return M


# --------------------------------------------------------- Jucys-Murphy


def jucys_murphy(lam: tuple[int, ...], k: int):
    """rho^lam(X_k) with X_k = sum_{i<k} (i k), 1-indexed:  X_k = sum over the
    transpositions swapping entries i and k for i = 1..k-1.

    In Young's orthogonal basis this is DIAGONAL, with the entry for tableau T
    equal to the content of the box holding k.
    """
    lam = tuple(x for x in lam if x > 0)
    n = sum(lam)
    f = dim_specht(lam)
    total = np.zeros((f, f))
    for i in range(1, k):
        # 0-indexed transposition swapping i-1 and k-1
        p = list(range(n))
        p[i - 1], p[k - 1] = p[k - 1], p[i - 1]
        total += young_orthogonal(lam, tuple(p))
    return total


def jucys_murphy_contents(lam, k):
    """The predicted diagonal of jucys_murphy(lam, k): content of the box
    holding the entry k, one value per standard tableau."""
    lam = tuple(x for x in lam if x > 0)
    out = []
    for T in syt_basis(lam):
        i, j = _positions(T)[k]
        out.append(j - i)
    return np.array(out, dtype=float)
