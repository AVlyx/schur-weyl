r"""Kronecker coefficients of the symmetric group.

Tensor two irreducible representations of $S_n$ and the result is no longer
irreducible. The Kronecker coefficient $g(\lambda, \mu, \nu)$ counts how many
copies of $S^\lambda$ sit inside $S^\mu \otimes S^\nu$. Since characters
multiply under tensor products, it is an inner product of characters: a finite
sum over conjugacy classes, with $z_\rho$ the centralizer size of the class
$\rho$,

$$
g(\lambda, \mu, \nu)
\;=\; \frac{1}{n!} \sum_{\sigma \in S_n} \chi^\lambda(\sigma)\, \chi^\mu(\sigma)\, \chi^\nu(\sigma)
\;=\; \sum_{\rho \vdash n} \frac{\chi^\lambda(\rho)\, \chi^\mu(\rho)\, \chi^\nu(\rho)}{z_\rho} .
$$

The same formula works for any number of partitions of the same $n$, so
`kronecker_coefficient` is variadic. With two arguments it reduces to character
orthogonality, $g(\lambda, \mu) = \delta_{\lambda\mu}$. No positive combinatorial
formula for these numbers is known in general — computing them from characters
is exactly what this does.

Examples:
    $S^{(2,1)} \otimes S^{(2,1)} = S^{(3)} \oplus S^{(2,1)} \oplus S^{(1,1,1)}$:

    >>> from schur_weyl import partitions, kronecker_coefficient
    >>> [kronecker_coefficient(lam, (2, 1), (2, 1)) for lam in partitions(3)]
    [1, 1, 1]
"""

from functools import lru_cache
from .symmetric_group import _z_for_permutation_conjugacy_class_size, permutation_conjugacy_classes
from .character import character
from math import prod


@lru_cache
def kronecker_coefficient(*partitions_: tuple[int, ...]) -> int:
    r"""Get the kronecker coefficient for multiple partitions.

    Args:
        *partitions_ (tuple[tuple[int, ...]]): The partitions

    Returns:
        int: The kronecker coefficient g(\mu, \nu, ...)

    Examples:
        >>> kronecker_coefficient((2, 1), (2, 1))
        1
        >>> kronecker_coefficient((2, 1), (3,))
        0
        >>> kronecker_coefficient((3,), (2, 1), (2, 1))
        1
        >>> kronecker_coefficient((3,), (2, 1), (1, 1, 1))
        0
    """
    if not partitions_:
        raise ValueError("Provide at least one partition")
    k: int = sum(partitions_[0])
    if not all(sum(partition) == k for partition in partitions_):
        raise ValueError("All partitions_ must belong to the same symmetric group Sn")
    ret = 0.0
    for cycle_type in permutation_conjugacy_classes(k):
        ret += prod(character(partition, cycle_type) for partition in partitions_) / _z_for_permutation_conjugacy_class_size(cycle_type)
    return round(ret)
