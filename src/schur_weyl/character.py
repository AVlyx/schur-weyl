"""Characters of the symmetric group via the Murnaghan-Nakayama rule."""

# Made by Claude

from functools import lru_cache
from .young_diagrams import partitions

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


def character(lam: tuple[int, ...], cycle_type: tuple[int, ...]) -> int:
    """chi^lam(cycle_type): the irreducible character of S_n indexed by the partition
    lam, evaluated on the conjugacy class of cycle type cycle_type.

    Examples:
        >>> character((3,), (1, 1, 1))       # trivial rep
        1
        >>> character((1, 1, 1), (2, 1))    # sign rep at a transposition
        -1
        >>> character((2, 1), (1, 1, 1))    # dimension of the standard rep
        2
        >>> character((2, 1), (3,))
        -1
        >>> character((2, 1), (2, 1))
        0
    """
    cycle_type = tuple(sorted((x for x in cycle_type if x > 0), reverse=True))
    if sum(lam) != sum(cycle_type):
        raise ValueError(f"size mismatch: |{lam}| != |{cycle_type}|")
    return _character(lam, cycle_type)


def character_table(n: int):
    """(rows, cols, table) with rows/cols the partitions of n in the order
    produced by partitions(n), and table[a][b] = chi^{rows[a]}(cols[b])."""
    ps = list(partitions(n))
    table = [[character(lam, mu) for mu in ps] for lam in ps]
    return ps, ps, table
