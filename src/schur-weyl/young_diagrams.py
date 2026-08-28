from typing import Generator, Optional


def partitions(k: int, max_part: Optional[int] = None, max_height: Optional[int] = None) -> Generator[tuple[int, ...]]:
    """Generate the partitions of k. All the non-increasing vectors of positive integers summing to k.

    Args:
        k (int): k
        max_part (Optional[int], optional): The maximum number of columns in any row of the Young tableau. Defaults to None.
        max_height (Optional[int], optional): The maximum number of rows. Defaults to None.

    Yields:
        Generator[tuple[int, ...]]: The partitions.

    Examples:
        >>> list(partitions(4))
        [(4,), (3, 1), (2, 2), (2, 1, 1), (1, 1, 1, 1)]

        >>> list(partitions(5, 2))
        [(2, 2, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)]

        >>> list(partitions(5, 2, 4))
        [(2, 2, 1), (2, 1, 1, 1)]
    """
    if max_part is None:
        max_part = k
    if max_height is None:
        max_height = k
    if k == 0:
        yield ()
        return
    if max_height == 0:
        return

    for first in range(min(k, max_part), 0, -1):
        for rest in partitions(k - first, first, max_height - 1):
            yield (first,) + rest


def partition_conjugate(lam: tuple[int, ...]) -> tuple[int, ...]:
    """Take the conjugate of a partition. Flip/Transpose a Young tableau accross its main diagonal

    Args:
        lam (tuple[int, ...]): The partition

    Returns:
        tuple[int, ...]: Its conjugate

    Examples:
        >>> partition_conjugate((4, 3, 1))
        (3, 2, 2, 1)
    """
    result: list[int] = [0] * lam[0]
    for p in lam:
        for j in range(p):
            result[j] += 1
    return tuple(result)


def hook_lenght(lam: tuple[int, ...], i: int, j: int) -> int:
    """The hook lenght at index ()

    Args:
        lam (tuple[int, ...]): The partition
        i (int): the row index (0-indexed)
        j (int): the column index (0-indexed)
    Returns:
        int: The hook lenght

    Examples:
        >>> hook_lenght((4, 3, 1), 0, 0)
        6
        >>> hook_lenght((4, 3, 3, 2), 2, 1)
        3
    """
    return (lam[i] - j) + sum(1 for l in range(i + 1, len(lam)) if lam[l] > j)  # (arm + 1) + leg


def cells(lam: tuple[int, ...]) -> Generator[tuple[int, int]]:
    """All the cells index of a partition

    Args:
        lam (tuple[int, ...]): The partition

    Yields:
        Generator[tuple[int, int]]: the indexes (0-indexed)

    Examples:
        >>> list(cells((3, 2, 1)))
        [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (2, 0)]
    """
    for i in range(len(lam)):
        for j in range(lam[i]):
            yield (i, j)


def majorizes(lam: tuple[int, ...], nu: tuple[int, ...]) -> bool:
    """Check if lam majorizes nu

    Args:
        lam (tuple[int, ...]): a partition
        nu (tuple[int, ...]): a partition of the same size as lambda. Otherwise the result is undefined

    Returns:
        bool: True if lam majorizes nu

    Examples:
        >>> majorizes((3, 3), (4, 1, 1))
        False
        >>> majorizes((4, 1, 1), (3, 3))
        False
        >>> majorizes((3, 1, 1, 1), (2, 2, 1, 1))
        True
        >>> majorizes((2, 1), (2, 1))
        True
    """
    if len(nu) < len(lam):
        return False
    c_lam, c_nu = 0, 0
    for i in range(len(lam)):
        c_lam += lam[i]
        c_nu += nu[i]
        if not c_lam >= c_nu:
            return False
    return True


def addable_corners(lam: tuple[int, ...]) -> Generator[tuple[int, int]]:
    """The cells that can be added to still obtain a valid partition

    Args:
        lam (tuple[int, ...]): The partition

    Yields:
        Generator[tuple[int, int]]: The indexes of the cells

    Examples:
        >>> list(addable_corners((2, 2, 1)))
        [(0, 2), (2, 1), (3, 0)]

    """
    if not lam:
        yield (0, 0)
        return
    yield 0, lam[0]
    for i in range(1, len(lam)):
        if lam[i] != lam[i - 1]:
            yield i, lam[i]
    yield len(lam), 0


def removable_corners(lam: tuple[int, ...]):
    """The cells that can be removed to still obtain a valid partition

    Args:
        lam (tuple[int, ...]): The partition

    Yields:
        Generator[tuple[int, int]]: The indexes of the cells

    Examples:
        >>> list(removable_corners((2, 2, 1)))
        [(2, 0), (1, 1)]
    """
    if not lam:
        return
    yield len(lam) - 1, lam[-1] - 1
    for i in range(len(lam) - 2, -1, -1):
        if lam[i] != lam[i + 1]:
            yield i, lam[i] - 1


def is_in_partition(lam: tuple[int, ...], i: int, j: int) -> bool:
    """Returns true if a cell is in the partition lam

    Args:
        lam (tuple[int, ...]): the partition
        i (int): the row index (0-indexed)
        j (int): the column index (0-indexed)

    Returns:
        bool: true if (i, j) in lam

    Examples:
        >>> is_in_partition((2, 1), 1, 1)
        False
        >>> is_in_partition((2, 1), 1, 0)
        True
    """
    return len(lam) > i and lam[i] > j
