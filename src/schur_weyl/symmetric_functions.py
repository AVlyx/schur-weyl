from collections import Counter

import numpy as np
from .tableaux import inverse_kostka_matrix


def schur_polynomial(lam: tuple[int, ...], xs: list[float]) -> float:
    """Get the schur polynomial of a partition evaluated on input xs
    Not very precise above len(lam) = 7

    Args:
        lam (tuple[int, ...]): The partition
        xs (list[float]): The input variables

    Returns:
        float: the value of the polynomial
    """
    return _jacobi_trudi(lam, xs)


def _jacobi_trudi(lam: tuple[int, ...], xs: list[float]) -> float:
    """Not very precise above len(lam) = 7. Could use an other det formula"""
    if not lam:
        return 1.0
    ell: int = len(lam)
    h_k: list[float] = h_0_to_k(ell + lam[0] - 1, xs)
    mat = np.array(
        [[h_k[lam[i] - i + j] if lam[i] - i + j >= 0 else 0.0 for j in range(ell)] for i in range(ell)],
        dtype=float,
    )
    return float(np.linalg.det(mat))


def h_0_to_k(k: int, xs: list[float]) -> list[float]:
    """Completely homogenous symmetric polynomial from h_0 to h_k

    Args:
        k (int): k
        xs (list[float]): The input values

    Returns:
        list[float]: h_0 to h_k
    """
    h_i: list[float] = [1.0] + [0.0] * k
    for xi in xs:
        for i in range(1, k + 1):
            h_i[i] += xi * h_i[i - 1]
    return h_i


def power_sum(mu: tuple[int, ...], xs: list[float]) -> float:
    """p_mu(x) = prod_i (sum_j x_j^{mu_i}).

    Examples:
        >>> power_sum((1,), [2, 3])        # p_1 = 2 + 3
        5
        >>> power_sum((2,), [2, 3])        # p_2 = 4 + 9
        13
        >>> power_sum((2, 1), [2, 3])     # p_2 * p_1 = 13 * 5
        65
        >>> power_sum((), [2, 3])
        1
    """
    result = 1
    for part in mu:
        result *= sum(xj**part for xj in xs)
    return result


def monomial_symmetric(mu: tuple[int, ...], xs: list[float]) -> float:
    """The monomial symmetric polynomial m_mu evaluated on input xs:
    the sum of x^alpha over the distinct rearrangements alpha of mu (padded with zeros to len(xs)).

    Computed variable by variable, keeping track of which parts of mu are still unassigned.

    Args:
        mu (tuple[int, ...]): The partition
        xs (list[float]): The input variables

    Returns:
        float: the value of the polynomial

    Examples:
        >>> monomial_symmetric((1, 1), [2, 3, 5])   # 2*3 + 2*5 + 3*5
        31
        >>> monomial_symmetric((2, 1), [2, 3])      # 4*3 + 9*2
        30
        >>> monomial_symmetric((1, 1, 1), [2, 3])   # more parts than variables
        0
        >>> monomial_symmetric((), [2, 3])
        1
    """
    counts = Counter(mu)
    parts = list(counts)
    # remaining multiplicity of each distinct part -> accumulated value
    state: dict[tuple[int, ...], float] = {tuple(counts[p] for p in parts): 1}
    for xi in xs:
        new_state: dict[tuple[int, ...], float] = {}
        for remaining, value in state.items():
            new_state[remaining] = new_state.get(remaining, 0) + value
            for j, part in enumerate(parts):
                if remaining[j]:
                    nxt = remaining[:j] + (remaining[j] - 1,) + remaining[j + 1 :]
                    new_state[nxt] = new_state.get(nxt, 0) + value * xi**part
        state = new_state
    return state.get((0,) * len(parts), 0)


def monomial_in_schur_basis(mu: tuple[int, ...]) -> list[tuple[int, tuple[int, ...]]]:
    """Expand the monomial symmetric function m_mu in the Schur basis.

    A row of inverse_kostka_matrix.

    Args:
        mu (tuple[int, ...]): The partition

    Returns:
        list[tuple[int, tuple[int, ...]]]: pairs (coefficient, nu) with m_mu = sum coefficient * s_nu

    Examples:
        >>> monomial_in_schur_basis((2,))
        [(1, (2,)), (-1, (1, 1))]
        >>> monomial_in_schur_basis((2, 1))
        [(1, (2, 1)), (-2, (1, 1, 1))]
        >>> monomial_in_schur_basis((3,))
        [(1, (3,)), (-1, (2, 1)), (1, (1, 1, 1))]
    """
    ps, table = inverse_kostka_matrix(sum(mu))
    row = table[ps.index(tuple(mu))]
    return [(c, nu) for c, nu in zip(row, ps) if c]
