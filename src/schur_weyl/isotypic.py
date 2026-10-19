r"""Isotypic projectors on $(\mathbb{C}^d)^{\otimes n}$.

The isotypic projector cuts $(\mathbb{C}^d)^{\otimes n}$ down to the $\lambda$
block. It is a class sum of the permutation action $U(\sigma)$, which shuffles
the tensor factors:

$$
P_\lambda \;=\; \frac{f^\lambda}{n!} \sum_{\sigma \in S_n} \chi^\lambda(\sigma)\, U(\sigma) .
$$

`isotypic_proj` builds this as a dense $d^n \times d^n$ matrix, which gets
expensive fast. When you only need to *apply* it, `apply_isotypic_proj` sums the
permuted copies of the tensor directly, so nothing larger than the input is
ever allocated. The $n$ axes you name are the tensor factors, in slot order;
every other axis is carried along untouched, which lets you project a batch of
tensors, or one leg of a larger network, in a single call.

Examples:
    >>> import numpy as np
    >>> from schur_weyl import partitions, isotypic_proj, apply_isotypic_proj
    >>> P = isotypic_proj((2, 1), d=2)
    >>> P.shape, bool(np.allclose(P @ P, P))
    ((8, 8), True)

    Axes 0 and 4 ride along, and the blocks resolve the identity:

    >>> T = np.random.default_rng(0).normal(size=(4, 2, 2, 2, 5))
    >>> apply_isotypic_proj(T, (2, 1), axes=(1, 2, 3)).shape
    (4, 2, 2, 2, 5)
    >>> bool(np.allclose(sum(apply_isotypic_proj(T, lam, axes=(1, 2, 3)) for lam in partitions(3)), T))
    True
"""

from __future__ import annotations

import numpy as np
from .symmetric_group import all_permutations_by_cycle_type, permutation_inverse
from .character import character
from .dimensions import dim_specht
from math import factorial


def _perm_index_map(sigma, d):
    n = len(sigma)
    idx = np.arange(d**n).reshape((d,) * n)
    return np.transpose(idx, axes=permutation_inverse(sigma)).reshape(-1)


def isotypic_proj(lam: tuple[int, ...], d: int) -> np.ndarray:
    """Return the isotypic projector acting on $(C^d)^{\\otimes k}$

    Args:
        lam (tuple[int, ...]): The partition k = sum(lam)
        d (int): The basis dimension

    Returns:
        np.ndarray: The projector
    """
    k = sum(lam)
    dimension = d**k
    rows = np.arange(dimension)

    proj = np.zeros((dimension, dimension), dtype=np.int64)

    for cycle, perms in all_permutations_by_cycle_type(k).items():
        char = character(lam, cycle)
        if char == 0:
            continue
        for perm in perms:
            proj[rows, _perm_index_map(perm, d)] += char

    return proj * (dim_specht(lam) / factorial(k))


def _axes_transpose(ndim: int, axes: tuple[int, ...], sigma: tuple[int, ...]) -> tuple[int, ...]:
    inv = permutation_inverse(sigma)
    perm = list(range(ndim))
    for i, axis in enumerate(axes):
        perm[axis] = axes[inv[i]]
    return tuple(perm)


def apply_isotypic_proj(tensor: np.ndarray, lam: tuple[int, ...], axes: tuple[int, ...] | None = None) -> np.ndarray:
    """Apply the isotypic projector $P_\\lambda$ along the chosen axes of a tensor.

    Args:
        tensor (np.ndarray): The tensor to project.
        lam (tuple[int, ...]): The partition
        axes (tuple[int, ...] | None): The k axes carrying the tensor factors. Defaults to the first k axes.

    Returns:
        np.ndarray: The projected tensor, of the same shape as tensor.

    Examples:
        >>> import numpy as np
        >>> v = np.arange(8.0).reshape(2, 2, 2)
        >>> p = apply_isotypic_proj(v, (2, 1))
        >>> bool(np.allclose(apply_isotypic_proj(p, (2, 1)), p))
        True
    """
    k = sum(lam)
    if not k:
        return tensor
    if axes is None:
        axes = tuple(range(k))
    axes = tuple(a % tensor.ndim for a in axes)
    if len(axes) != k:
        raise ValueError(f"expected {k} axes for lam={lam}, got {len(axes)}")
    if len(set(axes)) != k:
        raise ValueError(f"axes must be distinct, got {axes}")
    if len({tensor.shape[a] for a in axes}) != 1:
        raise ValueError(f"axes must have equal length, got {[tensor.shape[a] for a in axes]}")

    out = np.zeros(tensor.shape, dtype=np.result_type(tensor.dtype, np.float64))

    for cycle, perms in all_permutations_by_cycle_type(k).items():
        char = character(lam, cycle)
        if char == 0:
            continue
        acc = np.zeros_like(out)
        for perm in perms:
            acc += np.transpose(tensor, _axes_transpose(tensor.ndim, axes, perm))
        out += char * acc

    return out * (dim_specht(lam) / factorial(k))
