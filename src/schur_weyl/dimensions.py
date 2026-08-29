from math import factorial
from .young_diagrams import cells, hook_lenght


def dim_specht(lam: tuple[int, ...]) -> int:
    """The dimension of the Specht module S^lam.
    The number of SYT of shape lam.
    f^lam = n! / prod(hooks).

    Args:
        lam (tuple[int, ...]): The partition

    Returns:
        int: $f^\\lambda$

    Examples:
        >>> dim_specht((2, 1))
        2
        >>> dim_specht((3, 2, 1))
        16
    """
    n = sum(lam)
    hooks = 1
    for i, j in cells(lam):
        hooks *= hook_lenght(lam, i, j)
    return factorial(n) // hooks


def dim_weyl(lam: tuple[int, ...], d: int) -> int:
    """The dimension of the Weyl module V_lam for GL_d.
    The number of SSYT with entries from 1 to d

    Returns 0 if lam has more than d rows (a cell then has d + content = 0).

    Args:
        lam (tuple[int, ...]): The partition
        d (int): The alphabet size

    Returns:
        int: $dim V_lam^d$

    Examples:
        >>> dim_weyl((2, 1), 2)
        2
        >>> dim_weyl((2, 1), 3)
        8
        >>> dim_weyl((2,), 3)
        6
    """
    numerator = 1
    denominator = 1
    for i, j in cells(lam):
        numerator *= d + j - i
        denominator *= hook_lenght(lam, i, j)
    return numerator // denominator
