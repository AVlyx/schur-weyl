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
