from math import factorial
import pytest
from schur_weyl.dimensions import dim_specht, dim_weyl
from schur_weyl.young_diagrams import partitions


@pytest.mark.parametrize("n", range(10))
def test_sum_of_squared_specht_dims_is_n_factorial(n):
    assert sum(dim_specht(lam) ** 2 for lam in partitions(n)) == factorial(n)


@pytest.mark.parametrize("d", [2, 3, 4, 5])
@pytest.mark.parametrize("n", range(10))
def test_schur_weyl_dimension_count(n, d):
    assert sum(dim_specht(lam) * dim_weyl(lam, d) for lam in partitions(n, max_height=d)) == d**n
