import numpy as np


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
