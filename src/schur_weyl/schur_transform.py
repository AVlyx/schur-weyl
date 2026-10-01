"""The Schur transform: an explicit unitary realising Schur-Weyl duality,

    U : (C^d)^{otimes n} -> (+)_lam  V_lam^d (x) S^lam,

with the S_n factor in Young's orthogonal basis (`young_orthonormal`).

For every lam |- n with len(lam) <= d and every standard tableau T of shape lam
there is an isometry V_{lam,T} : C^{dim_weyl(lam, d)} -> (C^d)^{otimes n} with

    P(sigma) V_{lam,T} = sum_{T'} rho^lam(sigma)[T', T] V_{lam,T'}       (S_n acts on the tableau index)
    V_{lam,T}^dagger g^{otimes n} V_{lam,T'} = delta_{T,T'} q_lam(g)     (GL_d acts on the columns)

where P(sigma) is the permutation action of `isotypic._perm_index_map` (the tensor factor at
position i moves to position sigma(i)) and rho^lam = `young_orthogonal`.

Construction (deterministic, no eigensolver, everything real):
  1. T0 is the row-reading tableau of lam (1..lam_1 in the first row, and so on).  The joint
     eigenspace of the Jucys-Murphy elements X_j = sum_{i<j} P((i j)) with eigenvalues cont_T0(j)
     is projected onto with the product  prod_j prod_c (X_j - c) / (cont_T0(j) - c),  c running over
     the contents of the other addable corners of the shape of T0 restricted to 1..j-1.
  2. For each SSYT S of shape lam with entries in 1..d (in `semi_standard_young_tableau` order),
     w_S is the word with letter S(i, j) at position T0(i, j), i.e. the row reading word of S.
     The columns of V_{lam,T0} are the Gram-Schmidt orthonormalisation of P_T0 e_{w_S}.
     Sign convention: <e_{w_S}, column S of V_{lam,T0}> > 0.  Column S is a weight vector of
     weight content(S) under the diagonal torus.
  3. The other tableaux are reached by adjacent transpositions s_k:
        V_{s_k T} = (P(s_k) V_T - r V_T) / sqrt(1 - r^2),   r = 1 / axial_distance(T, k),
     which is exactly the column of Young's orthogonal form, so the bases are aligned.
"""

from __future__ import annotations

from collections import deque
from functools import lru_cache

import numpy as np

from .dimensions import dim_weyl
from .isotypic import _perm_index_map
from .tableaux import semi_standard_young_tableau
from .young_diagrams import addable_corners, partitions
from .young_orthonormal import _swap_entries, axial_distance, syt_basis


def _transposition(n: int, i: int, j: int) -> tuple[int, ...]:
    """The 0-indexed permutation swapping i and j."""
    p = list(range(n))
    p[i], p[j] = p[j], p[i]
    return tuple(p)


def _row_reading_tableau(lam: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    """The SYT filled with 1..n row by row."""
    out, start = [], 1
    for part in lam:
        out.append(tuple(range(start, start + part)))
        start += part
    return tuple(out)


def _jucys_murphy_projection(vectors: np.ndarray, T0, d: int) -> np.ndarray:
    """Project the columns of `vectors` onto the joint Jucys-Murphy eigenspace of T0."""
    n = sum(len(row) for row in T0)
    pos = {v: (i, j) for i, row in enumerate(T0) for j, v in enumerate(row)}
    shape: list[int] = []  # the shape of T0 restricted to the entries 1..k-1
    for k in range(1, n + 1):
        i, j = pos[k]
        target = j - i
        maps = [_perm_index_map(_transposition(n, a, k - 1), d) for a in range(k - 1)]
        for ci, cj in addable_corners(tuple(shape)):
            c = cj - ci
            if c == target:
                continue
            xv = sum((vectors[m] for m in maps), np.zeros_like(vectors))
            vectors = (xv - c * vectors) / (target - c)
        if i == len(shape):
            shape.append(0)
        shape[i] += 1
    return vectors


def _highest_tableau_block(lam: tuple[int, ...], d: int) -> np.ndarray:
    """V_{lam,T0}: the orthonormal basis of the T0 Jucys-Murphy eigenspace, see the module doc."""
    n = sum(lam)
    T0 = _row_reading_tableau(lam)
    positions = [v - 1 for row in T0 for v in row]
    ssyts = list(semi_standard_young_tableau(lam, d))
    words = np.zeros((d**n, len(ssyts)))
    strides = d ** np.arange(n - 1, -1, -1)
    for col, S in enumerate(ssyts):
        word = np.zeros(n, dtype=int)
        word[positions] = [letter - 1 for row in S for letter in row]
        words[int(word @ strides), col] = 1.0
    projected = _jucys_murphy_projection(words, T0, d)
    q, r = np.linalg.qr(projected)
    diag = np.diag(r)
    if np.abs(diag).min(initial=np.inf) < 1e-9:
        raise RuntimeError(f"projected SSYT words are not independent for lam={lam}, d={d}")
    return q * np.sign(diag)


@lru_cache(maxsize=None)
def schur_transform(d: int, n: int) -> dict[tuple[int, ...], np.ndarray]:
    """The Schur transform of (C^d)^{otimes n}, as one isometry per (lam, standard tableau).

    Args:
        d (int): The local dimension
        n (int): The number of tensor factors

    Returns:
        dict[tuple[int, ...], np.ndarray]: for each lam |- n with len(lam) <= d, a read-only
            array V of shape (dim_specht(lam), d**n, dim_weyl(lam, d)) with V[a] = V_{lam,T_a},
            T_a = syt_basis(lam)[a] the tableau indexing row/column a of `young_orthogonal`.
            Column s of every V[a] corresponds to semi_standard_young_tableau(lam, d)[s].

    Examples:
        >>> import numpy as np
        >>> V = schur_transform(2, 2)
        >>> {lam: block.shape for lam, block in V.items()}
        {(2,): (1, 4, 3), (1, 1): (1, 4, 1)}
        >>> np.round(V[(1, 1)][0, :, 0] * np.sqrt(2), 12) + 0.0   # the singlet
        array([ 0.,  1., -1.,  0.])
    """
    out: dict[tuple[int, ...], np.ndarray] = {}
    for lam in partitions(n):
        if len(lam) > d:
            continue
        basis = syt_basis(lam)
        index = {T: a for a, T in enumerate(basis)}
        blocks = np.zeros((len(basis), d**n, dim_weyl(lam, d)))
        T0 = _row_reading_tableau(lam)
        blocks[index[T0]] = _highest_tableau_block(lam, d)
        seen = {T0}
        queue = deque([T0])
        while queue:
            T = queue.popleft()
            for k in range(1, n):
                delta = axial_distance(T, k)
                if abs(delta) == 1:
                    continue
                S = _swap_entries(T, k)
                if S in seen:
                    continue
                r = 1.0 / delta
                VT = blocks[index[T]]
                blocks[index[S]] = (VT[_perm_index_map(_transposition(n, k - 1, k), d)] - r * VT) / np.sqrt(1.0 - r * r)
                seen.add(S)
                queue.append(S)
        blocks.flags.writeable = False
        out[lam] = blocks
    return out


def schur_transform_matrix(d: int, n: int) -> np.ndarray:
    """The full Schur transform as a real orthogonal d^n x d^n matrix U.

    The rows are grouped by lam in `partitions(n)` order; inside a block, row (s, a) =
    s * dim_specht(lam) + a is column s of V_{lam,T_a}.  Hence

        U P(sigma) U^T = (+)_lam  1_{dim_weyl(lam, d)} (x) young_orthogonal(lam, sigma)
        U g^{otimes n} U^T = (+)_lam  q_lam(g) (x) 1_{dim_specht(lam)}

    Examples:
        >>> import numpy as np
        >>> U = schur_transform_matrix(2, 3)
        >>> bool(np.allclose(U @ U.T, np.eye(8)))
        True
    """
    rows = [np.transpose(V, (2, 0, 1)).reshape(-1, d**n) for V in schur_transform(d, n).values()]
    return np.concatenate(rows, axis=0)
