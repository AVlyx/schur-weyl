"""Characters of the symmetric group via the Murnaghan-Nakayama rule."""

# Made by Claude

from functools import lru_cache
from typing import Generator, Optional
from math import factorial
from collections import Counter


def partitions(k: int, max_part: Optional[int] = None, max_height: Optional[int] = None) -> Generator[list[int], None, None]:
    if max_part is None:
        max_part = k
    if max_height is None:
        max_height = k
    if k == 0:
        yield []
        return
    if max_height == 0:
        return
    for first in range(min(k, max_part), 0, -1):
        for rest in partitions(k - first, first, max_height - 1):
            yield [first] + rest


def partition_conjugate(lam: list[int]) -> list[int]:
    if not lam:
        return []
    result = [0] * lam[0]
    for p in lam:
        for j in range(p):
            result[j] += 1
    return result


def cells(lam):
    for i in range(len(lam)):
        for j in range(lam[i]):
            yield (i, j)


def hook_lenght(lam: list[int], i: int, j: int) -> int:
    return (lam[i] - j) + sum(1 for l in range(i + 1, len(lam)) if lam[l] > j)


def dim_specht(lam: list[int]) -> int:
    n = sum(lam)
    hooks = 1
    for i, j in cells(lam):
        hooks *= hook_lenght(lam, i, j)
    return factorial(n) // hooks


def z_mu(mu) -> int:
    """The centraliser order: prod_i i^{m_i} m_i!."""
    z = 1
    for i, mi in Counter(mu).items():
        z *= (i**mi) * factorial(mi)
    return z


def conjugacy_class_size(mu) -> int:
    return factorial(sum(mu)) // z_mu(mu)


# ------------------------------------------------------------- border strips


def border_strips(lam: tuple[int, ...], k: int):
    """Yield (new_shape, height) for every border strip of size k removable
    from lam.

    Uses beta numbers (first-column hook lengths).  With L = len(lam), set
        beta_i = lam_i + (L - 1 - i),
    a strictly decreasing sequence of L distinct non-negative integers.
    Removing a border strip of size k is exactly: pick i with
    beta_i - k >= 0 and beta_i - k not already a beta value, then replace
    beta_i by beta_i - k.  The height of the strip is the number of beta
    values jumped over.
    """
    L = len(lam)
    if L == 0:
        return
    beta = [lam[i] + (L - 1 - i) for i in range(L)]
    bset = set(beta)
    for i in range(L):
        target = beta[i] - k
        if target < 0 or target in bset:
            continue
        height = sum(1 for b in beta if target < b < beta[i])
        new_beta = sorted(beta[:i] + [target] + beta[i + 1 :], reverse=True)
        new_lam = [new_beta[j] - (L - 1 - j) for j in range(L)]
        while new_lam and new_lam[-1] == 0:
            new_lam.pop()
        yield tuple(new_lam), height


# ------------------------------------------------------- Murnaghan-Nakayama


@lru_cache(maxsize=None)
def _character(lam: tuple[int, ...], mu: tuple[int, ...]) -> int:
    if not mu:
        return 1 if not lam else 0
    k = mu[0]
    rest = mu[1:]
    total = 0
    for new_lam, height in border_strips(lam, k):
        total += (-1) ** height * _character(new_lam, rest)
    return total


def character(lam, mu) -> int:
    """chi^lam(mu): the irreducible character of S_n indexed by the partition
    lam, evaluated on the conjugacy class of cycle type mu.

    Examples:
        >>> character([3], [1, 1, 1])       # trivial rep
        1
        >>> character([1, 1, 1], [2, 1])    # sign rep at a transposition
        -1
        >>> character([2, 1], [1, 1, 1])    # dimension of the standard rep
        2
        >>> character([2, 1], [3])
        -1
        >>> character([2, 1], [2, 1])
        0
    """
    lam = tuple(x for x in lam if x > 0)
    mu = tuple(sorted((x for x in mu if x > 0), reverse=True))
    if sum(lam) != sum(mu):
        raise ValueError(f"size mismatch: |{lam}| != |{mu}|")
    return _character(lam, mu)


def character_table(n: int):
    """(rows, cols, table) with rows/cols the partitions of n in the order
    produced by partitions(n), and table[a][b] = chi^{rows[a]}(cols[b])."""
    ps = [tuple(p) for p in partitions(n)]
    table = [[character(lam, mu) for mu in ps] for lam in ps]
    return ps, ps, table
