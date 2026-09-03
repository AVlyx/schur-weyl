from itertools import permutations
from typing import Any, Generator
from collections import Counter, defaultdict
from math import factorial
from .young_diagrams import partitions


def permutation_identity(n) -> tuple[int, ...]:
    """The identity permutation on {0, ..., n-1}.

    Examples:
        >>> permutation_identity(4)
        (0, 1, 2, 3)
    """
    return tuple(range(n))


def permutation_compose(pi: tuple[int, ...], sigma: tuple[int, ...]) -> tuple[int, ...]:
    """Multiply two permutation (pi * sigma)(i) = pi(sigma(i))

    Examples:
        >>> permutation_compose((1, 0, 2), (0, 2, 1))
        (1, 2, 0)
        >>> permutation_compose((0, 1, 2), (2, 0, 1))
        (2, 0, 1)
    """
    return tuple(pi[sigma[i]] for i in range(len(pi)))


def permutation_inverse(sigma: tuple[int, ...]) -> tuple[int, ...]:
    """The inverse permutation.

    Examples:
        >>> permutation_inverse((1, 2, 0))
        (2, 0, 1)
        >>> permutation_inverse((0, 1, 2))
        (0, 1, 2)
    """
    result = [0] * len(sigma)
    for i, s in enumerate(sigma):
        result[s] = i
    return tuple(result)


def all_permutations(n: int) -> Generator[tuple[int, ...], None, None]:
    """All permutations of Sn

    Args:
        n (int): n

    Yields:
        Generator[tuple[int, ...], None, None]: the permutations
    """
    for p in permutations(range(n)):
        yield p


def all_permutations_by_cycle_type(n: int) -> dict[tuple[int, ...], list[tuple[int, ...]]]:
    """All permutations of Sn assigned to their cycle type

    Args:
        n (int): n

    Yields:
        dict[tuple[int, ...], list[tuple[int, ...]]]: the cycle type and its permutations
    """
    cycle_and_perm: dict[tuple[int, ...], list[tuple[int, ...]]] = defaultdict(list)
    for p in permutations(range(n)):
        cycle_and_perm[permutation_cycle_type(p)].append(p)
    return cycle_and_perm


def permutation_cycle_type(perm: tuple[int, ...]) -> tuple[int, ...]:
    """Give the cycle type of a permutation. The disjoint sets of perms

    Args:
        perm (tuple[int, ...]): The permutation

    Returns:
        tuple[int, ...]: Cycle type

    Examples:
        >>> permutation_cycle_type((1, 2, 0, 4, 3))
        (3, 2)
        >>> permutation_cycle_type((0, 1, 2))
        (1, 1, 1)
        >>> permutation_cycle_type((1, 0, 2))
        (2, 1)
    """

    seen: list[bool] = [False] * len(perm)
    cycles: list[int] = []
    for i in range(len(perm)):
        if seen[i]:
            continue
        seen[i] = True

        count = 1
        j = perm[i]
        while not seen[j]:
            seen[j] = True
            j = perm[j]
            count += 1
        cycles.append(count)
    return tuple(sorted(cycles, reverse=True))


def permutation_conjugacy_classes(k: int) -> Generator[tuple[int, ...], None, None]:
    """All the permutation conjugacy class represented by their cycle type.
    Simply the partitions of k.

    Args:
        k (int): The size of the Symmetric group

    Yields:
        Generator[tuple[int, ...], None, None]: _description_
    """
    return partitions(k)


def _z_for_permutation_conjugacy_class_size(cycle_type: tuple[int, ...]) -> int:
    counts = Counter(cycle_type)
    z: int = 1
    for cycle_len, cycle_count in counts.items():
        z *= cycle_len**cycle_count * factorial(cycle_count)
    return z


def permutation_conjugacy_class_size(cycle_type: tuple[int, ...]) -> int:
    """The number of permutations of the same cycle type

    Args:
        cycle_type (tuple[int, ...]): The cycle type

    Returns:
        int: The number of permutations of the same cycle type

    Examples:
        >>> permutation_conjugacy_class_size((1,1,1,1))
        1
        >>> permutation_conjugacy_class_size((2,1,1))
        6
    """

    return factorial(sum(cycle_type)) // _z_for_permutation_conjugacy_class_size(cycle_type)
