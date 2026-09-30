from functools import lru_cache
from typing import Generator

import numpy as np
from .young_diagrams import partitions, removable_corners


def standard_young_tableaux(lam: tuple[int, ...]) -> Generator[tuple[tuple[int, ...], ...], None, None]:
    """Generates all SYT of a given shape

    Args:
        lam (tuple[int, ...]): The partition

    Yields:
        Generator[tuple[tuple[int, ...], ...], None, None]: a SYT

    Examples:
        >>> list(standard_young_tableaux(()))
        [()]
        >>> sorted(standard_young_tableaux((2, 1)))
        [((1, 2), (3,)), ((1, 3), (2,))]
        >>> len(list(standard_young_tableaux((3, 2, 1))))
        16
    """
    if not lam:
        yield ()
        return

    ret: list[list[int]] = [[-1] * li for li in lam]
    ret[0][0] = 1

    def syt_rec(small_lam: list[int], d: int) -> Generator[tuple[tuple[int, ...], ...], None, None]:
        if len(small_lam) == 1 and small_lam[0] == 1:
            yield tuple(tuple(row) for row in ret)
            return
        for i, j in removable_corners(tuple(small_lam)):
            if j == 0:
                small_lam.pop(i)
                ret[i][j] = d
                yield from syt_rec(small_lam, d - 1)
                small_lam.append(1)
            else:
                small_lam[i] -= 1
                ret[i][j] = d
                yield from syt_rec(small_lam, d - 1)
                small_lam[i] += 1

    yield from syt_rec(list(lam), sum(lam))


def semi_standard_young_tableau(lam: tuple[int, ...], d: int) -> Generator[tuple[tuple[int, ...], ...], None, None]:
    """Generate all SSYT of a given shape and alphabet

    Args:
        lam (tuple[int, ...]): The partition
        d (int): the alphabet size

    Yields:
        Generator[tuple[tuple[int, ...], ...], None, None]: The SSYT

    Examples:
        >>> list(semi_standard_young_tableau((2, 1), 3))
        [((1, 1), (2,)), ((1, 2), (2,)), ((1, 3), (2,)), ((1, 1), (3,)), ((1, 2), (3,)), ((1, 3), (3,)), ((2, 2), (3,)), ((2, 3), (3,))]
    """
    if not lam:
        yield ()
        return
    ret: list[list[int]] = [[-1] * li for li in lam]
    maxv: list[list[int]] = [[d + 1] * li for li in lam]

    def fill(i: int, j: int, minv: int) -> Generator[tuple[tuple[int, ...], ...], None, None]:
        if j >= lam[i]:
            if i == 0:
                yield tuple(tuple(row) for row in ret)
                return
            yield from fill(i - 1, 0, i)
            return

        for cell_val in range(minv, maxv[i][j]):
            ret[i][j] = cell_val
            if i != 0:
                maxv[i - 1][j] = cell_val
            yield from fill(i, j + 1, cell_val)

    yield from fill(len(lam) - 1, 0, len(lam))


def semi_standard_young_tableau_by_content(lam: tuple[int, ...], content: tuple[int, ...]) -> Generator[tuple[tuple[int, ...], ...], None, None]:
    """Generate all SSYT with given content

    Args:
        lam (tuple[int, ...]): The partition
        content (tuple[int, ...]): The content (Should satisfy sum(lam) == sum(content))

    Yields:
        Generator[tuple[tuple[int, ...], ...], None, None]: a SSYT

    Examples:
        >>> list(semi_standard_young_tableau_by_content((2, 1), (2, 1)))
        [((1, 1), (2,))]
        >>> list(semi_standard_young_tableau_by_content((2, 1), (1, 1, 1)))
        [((1, 3), (2,)), ((1, 2), (3,))]
        >>> list(semi_standard_young_tableau_by_content((), ()))
        [()]
    """
    if not sum(lam) == sum(content):
        return
    if not lam:
        yield ()
        return

    d: int = len(content)
    contentcp: list[int] = list(content)
    ret: list[list[int]] = [[-1] * li for li in lam]
    maxv: list[list[int]] = [[d + 1] * li for li in lam]

    def fill(i: int, j: int, minv: int) -> Generator[tuple[tuple[int, ...], ...], None, None]:
        if j >= lam[i]:
            if i == 0:
                yield tuple(tuple(row) for row in ret)
                return
            yield from fill(i - 1, 0, i)
            return

        for cell_val in range(minv, maxv[i][j]):
            if contentcp[cell_val - 1] <= 0:
                continue
            ret[i][j] = cell_val
            if i != 0:
                maxv[i - 1][j] = cell_val
            contentcp[cell_val - 1] -= 1
            yield from fill(i, j + 1, cell_val)
            contentcp[cell_val - 1] += 1

    yield from fill(len(lam) - 1, 0, len(lam))


@lru_cache(maxsize=None)
def kostka(lam: tuple[int, ...], mu: tuple[int, ...]) -> int:
    """Compute the kostka number: the number of SSYT of shape lam and content mu.

    Recursion on horizontal strips: the cells filled with the largest letter form a
    horizontal strip of size mu[-1], so K_{lam, mu} = sum_nu K_{nu, mu[:-1]}.

    Args:
        lam (tuple[int, ...]): The partition
        mu (tuple[int, ...]): content, any composition (zeros allowed)

    Returns:
        int: $K_{\\lambda \\mu}$

    Examples:
        >>> kostka((2, 1), (1, 1, 1))
        2
        >>> kostka((2, 1), (0, 2, 1))
        1
        >>> from math import factorial
        >>> kostka((10, 1, 1, 1, 1, 1, 1, 1), (1,) * 17) == factorial(17) // (17 * factorial(9) * factorial(7))
        True
    """
    lam, mu = tuple(lam), tuple(mu)
    while lam and lam[-1] == 0:
        lam = lam[:-1]
    if sum(lam) != sum(mu):
        return 0
    if not mu:
        return 1
    if len(lam) > len(mu):  # the first column needs distinct letters
        return 0
    return sum(kostka(nu, mu[:-1]) for nu in _remove_horizontal_strips(lam, mu[-1]))


def _remove_horizontal_strips(lam: tuple[int, ...], k: int) -> Generator[tuple[int, ...], None, None]:
    """All nu with lam / nu a horizontal strip of size k, i.e. lam[i+1] <= nu[i] <= lam[i]."""
    nu: list[int] = [0] * len(lam)

    def fill(i: int, remaining: int) -> Generator[tuple[int, ...], None, None]:
        if i == len(lam):
            if remaining == 0:
                yield tuple(part for part in nu if part)
            return
        floor = lam[i + 1] if i + 1 < len(lam) else 0
        for removed in range(min(remaining, lam[i] - floor) + 1):
            nu[i] = lam[i] - removed
            yield from fill(i + 1, remaining - removed)

    yield from fill(0, k)


@lru_cache(maxsize=None)
def kostka_matrix(n: int) -> tuple[tuple[tuple[int, ...], ...], tuple[tuple[int, ...], ...]]:
    """All the kostka numbers K_{lam, mu} for lam, mu partitions of n.

    Pieri rule: h_mu = sum_lam K_{lam, mu} s_lam, so a column is built by adding horizontal
    strips of sizes mu[0], mu[1], ... to the empty shape. partitions(n) keeps partitions with a
    common prefix together, so the columns of the current prefixes are kept on a stack.

    Args:
        n (int): The size of the partitions

    Returns:
        tuple[tuple[tuple[int, ...], ...], tuple[tuple[int, ...], ...]]: (ps, table) with ps the partitions
            of n in the order produced by partitions(n), and table[a][b] = K_{ps[a], ps[b]}

    Examples:
        >>> kostka_matrix(3)
        (((3,), (2, 1), (1, 1, 1)), ((1, 1, 1), (0, 1, 2), (0, 0, 1)))
    """
    ps = tuple(partitions(n))
    index = {lam: a for a, lam in enumerate(ps)}
    table = [[0] * len(ps) for _ in ps]
    strips: dict[tuple[tuple[int, ...], int], list[tuple[int, ...]]] = {}
    stack: list[dict[tuple[int, ...], int]] = [{(): 1}]  # stack[i] is the column of mu[:i]
    prev: tuple[int, ...] = ()
    for b, mu in enumerate(ps):
        common = 0
        while common < min(len(prev), len(mu)) and prev[common] == mu[common]:
            common += 1
        del stack[common + 1 :]
        for k in mu[common:]:
            column: dict[tuple[int, ...], int] = {}
            for nu, count in stack[-1].items():
                if (nu, k) not in strips:
                    strips[nu, k] = list(_add_horizontal_strips(nu, k))
                for lam in strips[nu, k]:
                    column[lam] = column.get(lam, 0) + count
            stack.append(column)
        for lam, count in stack[-1].items():
            table[index[lam]][b] = count
        prev = mu
    return ps, tuple(tuple(row) for row in table)


@lru_cache(maxsize=None)
def inverse_kostka_matrix(n: int) -> tuple[tuple[tuple[int, ...], ...], tuple[tuple[int, ...], ...]]:
    """The inverse of kostka_matrix(n), which changes the Schur basis into the monomial basis:
    m_mu = sum_nu table[mu][nu] s_nu.

    kostka_matrix(n) is upper unitriangular, so the inverse is solved column by column with numpy,
    in int64 when that cannot overflow and with python ints otherwise.

    Args:
        n (int): The size of the partitions

    Returns:
        tuple[tuple[tuple[int, ...], ...], tuple[tuple[int, ...], ...]]: (ps, table) with ps the partitions
            of n in the order produced by partitions(n), and table[a][b] the entry of the inverse

    Examples:
        >>> inverse_kostka_matrix(3)
        (((3,), (2, 1), (1, 1, 1)), ((1, -1, 1), (0, 1, -2), (0, 0, 1)))
    """
    ps, table = kostka_matrix(n)
    inverse = _unitriangular_inverse(table, np.int64)
    if inverse is None:
        inverse = _unitriangular_inverse(table, object)
    assert inverse is not None
    return ps, tuple(tuple(row) for row in inverse.tolist())


def _unitriangular_inverse(table: tuple[tuple[int, ...], ...], dtype) -> np.ndarray | None:
    """Inverse of an upper unitriangular integer matrix, or None if int64 could overflow."""
    try:
        a = np.array(table, dtype=dtype).reshape(len(table), len(table))
    except OverflowError:
        return None
    x = np.eye(len(table), dtype=int).astype(dtype)
    exact = dtype is object
    abs_a = np.abs(a).astype(float)
    abs_x = x.astype(float)
    for c in range(1, len(table)):
        # x a = I column c: x[:c, c] = -x[:c, :c] a[:c, c], as x[c, c] = 1 and x is upper triangular
        if not exact and (abs_x[:c, :c] @ abs_a[:c, c]).max() >= 2.0**62:
            return None
        x[:c, c] = -(x[:c, :c] @ a[:c, c])
        if not exact:
            abs_x[:c, c] = np.abs(x[:c, c])
    return x


def _add_horizontal_strips(nu: tuple[int, ...], k: int) -> Generator[tuple[int, ...], None, None]:
    """All lam with lam / nu a horizontal strip of size k, i.e. nu[i] <= lam[i] <= nu[i-1]."""
    padded = nu + (0,)
    lam: list[int] = [0] * len(padded)

    def fill(i: int, remaining: int) -> Generator[tuple[int, ...], None, None]:
        if i == len(padded):
            if remaining == 0:
                yield tuple(part for part in lam if part)
            return
        ceiling = padded[i - 1] - padded[i] if i > 0 else remaining
        for added in range(min(remaining, ceiling) + 1):
            lam[i] = padded[i] + added
            yield from fill(i + 1, remaining - added)

    yield from fill(0, k)


def reverse_reading_word(tableau: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    """The reverse reading word of a tableau: right-to-left within each row,
    bottom row to top row.

    Examples:
        >>> reverse_reading_word(((1, 1, 3), (2, 3), (4,)))
        (4, 3, 2, 3, 1, 1)
        >>> reverse_reading_word(((1, 2), (3,)))
        (3, 2, 1)
        >>> reverse_reading_word(())
        ()
    """
    return tuple(v for row in reversed(tableau) for v in reversed(row))
