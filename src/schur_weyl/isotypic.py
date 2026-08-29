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
