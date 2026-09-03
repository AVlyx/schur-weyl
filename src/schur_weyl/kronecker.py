from .symmetric_group import _z_for_permutation_conjugacy_class_size, permutation_conjugacy_classes
from .character import character
from math import prod


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
