import pytest

from schur_weyl.young_diagrams import partition_conjugate, partitions

PARTITION_COUNTS = [1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42]  # n = 0..10


@pytest.mark.parametrize("n, count", list(enumerate(PARTITION_COUNTS)))
def test_partition_count(n, count):
    assert len(list(partitions(n))) == count


@pytest.mark.parametrize("n", range(11))
def test_conjugate_is_an_involution(n):
    for lam in partitions(n):
        assert partition_conjugate(partition_conjugate(lam)) == lam
