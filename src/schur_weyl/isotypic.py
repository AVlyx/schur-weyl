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
