from typing import Generator
from young_diagrams import removable_corners


def standard_young_tableaux(lam: list[int]) -> Generator[list[list[int]]]:
    """Generates all SYT of a given shape

    Args:
        lam (list[int]): The partition

    Yields:
        Generator[list[list[int]]]: a SYT

    Examples:
        >>> list(standard_young_tableaux([]))
        [[]]
        >>> sorted(standard_young_tableaux([2, 1]))
        [[[1, 2], [3]], [[1, 3], [2]]]
        >>> len(list(standard_young_tableaux([3, 2, 1])))
        16
    """
    if not lam:
        yield []
        return

    ret: list[list[int]] = [[-1] * li for li in lam]
    ret[0][0] = 1

    def syt_rec(small_lam: list[int], d: int) -> Generator[list[list[int]]]:
        if len(small_lam) == 1 and small_lam[0] == 1:
            yield [row.copy() for row in ret]
            return
        for i, j in removable_corners(small_lam):
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

    yield from syt_rec(lam, sum(lam))


def semi_standard_young_tableau(lam: list[int], d: int) -> Generator[list[list[int]]]:
    """Generate all SSYT of a given shape and alphabet

    Args:
        lam (list[int]): The partition
        d (int): the alphabet size

    Yields:
        Generator[list[list[int]]]: The SSYT

    Examples:
        >>> list(semi_standard_young_tableau([2, 1], 3))
        [[[1, 1], [2]], [[1, 2], [2]], [[1, 3], [2]], [[1, 1], [3]], [[1, 2], [3]], [[1, 3], [3]], [[2, 2], [3]], [[2, 3], [3]]]
    """
    if not lam:
        yield []
        return
    ret: list[list[int]] = [[-1] * li for li in lam]
    maxv: list[list[int]] = [[d + 1] * li for li in lam]

    def fill(i: int, j: int, minv: int) -> Generator[list[list[int]]]:
        if j >= lam[i]:
            if i == 0:
                yield [row.copy() for row in ret]
                return
            yield from fill(i - 1, 0, i)
            return

        for cell_val in range(minv, maxv[i][j]):
            ret[i][j] = cell_val
            if i != 0:
                maxv[i - 1][j] = cell_val
            yield from fill(i, j + 1, cell_val)

    yield from fill(len(lam) - 1, 0, len(lam))


def semi_standard_young_tableau_by_content(lam: list[int], content: list[int]) -> Generator[list[list[int]]]:
    """Generate all SSYT with given content

    Args:
        lam (list[int]): The partition
        content (list[int]): The content (Should satisfy sum(lam) == sum(content))

    Yields:
        Generator[list[list[int]]]: a SSYT

    Examples:
        >>> list(semi_standard_young_tableau_by_content([2, 1], [2, 1]))
        [[[1, 1], [2]]]
        >>> list(semi_standard_young_tableau_by_content([2, 1], [1, 1, 1]))
        [[[1, 3], [2]], [[1, 2], [3]]]
        >>> list(semi_standard_young_tableau_by_content([], []))
        [[]]
    """
    if not sum(lam) == sum(content):
        return
    if not lam:
        yield []
        return

    d: int = len(content)
    contentcp: list[int] = content.copy()
    ret: list[list[int]] = [[-1] * li for li in lam]
    maxv: list[list[int]] = [[d + 1] * li for li in lam]

    def fill(i: int, j: int, minv: int) -> Generator[list[list[int]]]:
        if j >= lam[i]:
            if i == 0:
                yield [row.copy() for row in ret]
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


# TODO implement with better algorithm (horizontal strips)
def kostka(lam: list[int], mu: list[int]) -> int:
    """Compute the kotska number

    Args:
        lam (list[int]): The partition
        mu (list[int]): content

    Returns:
        int: $K_{\\lambda \\mu}$

    Examples:
        >>> kostka([2, 1], [1, 1, 1])
        2
        >>> from math import factorial
        >>> kostka([10,1,1,1,1,1,1,1], [1] * 17) == factorial(17) / (17 * factorial(9) * factorial(7))
        True
    """
    return len(list(semi_standard_young_tableau_by_content(lam, mu)))


def reverse_reading_word(tableau: list[list[int]]) -> list[int]:
    """The reverse reading word of a tableau: right-to-left within each row,
    bottom row to top row.

    Examples:
        >>> reading_word([[1, 1, 3], [2, 3], [4]])
        [4, 3, 2, 3, 1, 1]
        >>> reading_word([[1, 2], [3]])
        [3, 2, 1]
        >>> reading_word([])
        []
    """
    return [v for row in reversed(tableau) for v in reversed(row)]
