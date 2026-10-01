from functools import reduce

import numpy as np
import pytest

from schur_weyl.dimensions import dim_specht, dim_weyl
from schur_weyl.isotypic import _perm_index_map, isotypic_proj
from schur_weyl.schur_transform import schur_transform, schur_transform_matrix
from schur_weyl.symmetric_functions import schur_polynomial
from schur_weyl.symmetric_group import permutation_compose
from schur_weyl.tableaux import semi_standard_young_tableau
from schur_weyl.young_diagrams import partitions
from schur_weyl.young_orthonormal import jucys_murphy_contents, syt_basis, young_orthogonal

CASES = [(1, 3), (2, 2), (2, 3), (2, 4), (2, 5), (3, 3), (3, 4), (4, 3), (4, 4), (2, 6), (3, 5)]


def perm_matrix(sigma, d):
    """P(sigma) on (C^d)^{otimes n}, as in tests/test_isotypic.py."""
    return np.eye(d ** len(sigma))[_perm_index_map(sigma, d)]


def transposition(n, i, j):
    p = list(range(n))
    p[i], p[j] = p[j], p[i]
    return tuple(p)


def random_perm(rng, n):
    return tuple(int(x) for x in rng.permutation(n))


def tensor_power(g, n):
    return reduce(np.kron, [g] * n, np.eye(1))


def irrep_block(V, g, n):
    """q_lam(g) read off from the first tableau."""
    return V[0].conj().T @ tensor_power(g, n) @ V[0]


# 1. shapes and dimensions


@pytest.mark.parametrize("d, n", CASES)
def test_shapes_and_dimensions(d, n):
    V = schur_transform(d, n)
    assert list(V) == [lam for lam in partitions(n) if len(lam) <= d]
    for lam, block in V.items():
        assert block.shape == (dim_specht(lam), d**n, dim_weyl(lam, d))
        assert len(syt_basis(lam)) == block.shape[0]
        assert block.shape[2] > 0
    assert sum(dim_specht(lam) * dim_weyl(lam, d) for lam in V) == d**n


# 2. unitarity


@pytest.mark.parametrize("d, n", CASES)
def test_unitarity(d, n):
    columns = [(lam, a, block[a]) for lam, block in schur_transform(d, n).items() for a in range(block.shape[0])]
    total = np.zeros((d**n, d**n))
    for lam, a, Va in columns:
        for mu, b, Wb in columns:
            expected = np.eye(Va.shape[1]) if (lam, a) == (mu, b) else 0
            assert np.allclose(Va.T @ Wb, expected)
        total += Va @ Va.T
    assert np.allclose(total, np.eye(d**n))
    U = schur_transform_matrix(d, n)
    assert np.allclose(U @ U.T, np.eye(d**n))


# 3. agreement with the projectors


@pytest.mark.parametrize("d, n", CASES)
def test_isotypic_projectors(d, n):
    for lam, block in schur_transform(d, n).items():
        assert np.allclose(sum(Va @ Va.T for Va in block), isotypic_proj(lam, d))


@pytest.mark.parametrize("d, n", [(2, 4), (3, 3), (4, 4)])
def test_symmetric_and_antisymmetric(d, n):
    V = schur_transform(d, n)
    sym = V[(n,)][0]
    for sigma in [transposition(n, 0, 1), tuple(range(1, n)) + (0,)]:
        assert np.allclose(perm_matrix(sigma, d) @ sym, sym)
    if n <= d:
        anti = V[(1,) * n][0]
        assert anti.shape[1] == dim_weyl((1,) * n, d)
        assert np.allclose(perm_matrix(transposition(n, 0, 1), d) @ anti, -anti)
        assert np.allclose(anti @ anti.T, isotypic_proj((1,) * n, d))


# 4. permutation action and aligned bases


@pytest.mark.parametrize("d, n", CASES)
def test_adjacent_transpositions(d, n):
    for lam, block in schur_transform(d, n).items():
        for k in range(n - 1):
            s = transposition(n, k, k + 1)
            rho = young_orthogonal(lam, s)
            P = perm_matrix(s, d)
            for a in range(block.shape[0]):
                assert np.allclose(P @ block[a], np.einsum("b,bij->ij", rho[:, a], block))


@pytest.mark.parametrize("d, n", CASES)
def test_random_permutations_and_block_form(d, n):
    rng = np.random.default_rng(n * 10 + d)
    U = schur_transform_matrix(d, n)
    V = schur_transform(d, n)
    sigmas = [random_perm(rng, n) for _ in range(5)]
    if n >= 3:
        sigmas += [(1, 2, 0) + tuple(range(3, n)), (2, 0, 1) + tuple(range(3, n))]  # both 3-cycles
    for sigma in sigmas:
        P = perm_matrix(sigma, d)
        for lam, block in V.items():
            rho = young_orthogonal(lam, sigma)
            for a in range(block.shape[0]):
                assert np.allclose(P @ block[a], np.einsum("b,bij->ij", rho[:, a], block))
        expected = [np.kron(np.eye(dim_weyl(lam, d)), young_orthogonal(lam, sigma)) for lam in V]
        assert np.allclose(U @ P @ U.T, _block_diag(expected))


def test_conventions_compose():
    # P and rho must both be homomorphisms with the same convention
    n, d = 4, 2
    rng = np.random.default_rng(1)
    for _ in range(5):
        s, t = random_perm(rng, n), random_perm(rng, n)
        st = permutation_compose(s, t)
        assert np.allclose(perm_matrix(s, d) @ perm_matrix(t, d), perm_matrix(st, d))
        for lam in partitions(n):
            assert np.allclose(young_orthogonal(lam, s) @ young_orthogonal(lam, t), young_orthogonal(lam, st))


@pytest.mark.parametrize("d, n", CASES)
def test_jucys_murphy_contents(d, n):
    V = schur_transform(d, n)
    for j in range(1, n + 1):
        X = sum((perm_matrix(transposition(n, i, j - 1), d) for i in range(j - 1)), np.zeros((d**n, d**n)))
        for lam, block in V.items():
            contents = jucys_murphy_contents(lam, j)
            for a in range(block.shape[0]):
                assert np.allclose(block[a].T @ X @ block[a], contents[a] * np.eye(block.shape[2]))


# 5. action of GL_d


@pytest.mark.parametrize("d, n", CASES)
def test_general_linear_action(d, n):
    rng = np.random.default_rng(d * 100 + n)
    g = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    h = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    gn, hn, ghn = tensor_power(g, n), tensor_power(h, n), tensor_power(g @ h, n)
    eig = np.linalg.eigvals(g)
    for lam, block in schur_transform(d, n).items():
        q = block[0].T @ gn @ block[0]
        for a in range(block.shape[0]):
            for b in range(block.shape[0]):
                assert np.allclose(block[a].T @ gn @ block[b], q if a == b else 0)
        assert np.isclose(np.trace(q), _schur_complex(lam, eig))
        qh = block[0].T @ hn @ block[0]
        assert np.allclose(block[0].T @ ghn @ block[0], q @ qh)


def _schur_complex(lam, xs):
    # schur_polynomial works on floats; compare a complex evaluation via the real and imaginary
    # parts of the bialternant formula instead
    m = len(xs)
    lam = tuple(lam) + (0,) * (m - len(lam))
    num = np.linalg.det(np.array([[x ** (lam[j] + m - 1 - j) for j in range(m)] for x in xs]))
    den = np.linalg.det(np.array([[x ** (m - 1 - j) for j in range(m)] for x in xs]))
    return num / den


@pytest.mark.parametrize("d, n", [(2, 4), (3, 3), (3, 4)])
def test_characters_against_schur_polynomial(d, n):
    rng = np.random.default_rng(7)
    a = rng.normal(size=(d, d))
    g = a @ a.T + np.eye(d)  # real symmetric, real spectrum
    gn = tensor_power(g, n)
    for lam, block in schur_transform(d, n).items():
        trace = np.trace(block[0].T @ gn @ block[0])
        assert np.isclose(trace, schur_polynomial(lam, list(np.linalg.eigvalsh(g))))


# 6. realness, determinism, sign convention


@pytest.mark.parametrize("d, n", CASES)
def test_real_and_deterministic(d, n):
    V = schur_transform(d, n)
    schur_transform.cache_clear()
    W = schur_transform(d, n)
    for lam in V:
        assert V[lam].dtype == np.float64
        assert np.array_equal(V[lam], W[lam])
        assert not V[lam].flags.writeable
    rng = np.random.default_rng(3)
    g = rng.normal(size=(d, d))
    for block in V.values():
        assert np.isrealobj(irrep_block(block, g, n))


@pytest.mark.parametrize("d, n", CASES)
def test_sign_convention(d, n):
    # <e_{w_S}, column S of V_{lam,T0}> > 0, w_S the row reading word of S, T0 = syt_basis(lam)[0]
    strides = d ** np.arange(n - 1, -1, -1)
    for lam, block in schur_transform(d, n).items():
        assert syt_basis(lam)[0] == _row_reading(lam)
        for s, S in enumerate(semi_standard_young_tableau(lam, d)):
            word = np.array([letter - 1 for row in S for letter in row], dtype=int)
            assert block[0][int(word @ strides) if n else 0, s] > 1e-9


def _row_reading(lam):
    out, start = [], 1
    for part in lam:
        out.append(tuple(range(start, start + part)))
        start += part
    return tuple(out)


# 7. edge cases


def test_edge_cases():
    assert {lam: b.tolist() for lam, b in schur_transform(3, 0).items()} == {(): [[[1.0]]]}
    V = schur_transform(3, 1)
    assert list(V) == [(1,)] and np.allclose(V[(1,)][0], np.eye(3))
    V = schur_transform(1, 5)
    assert list(V) == [(5,)] and np.allclose(V[(5,)], np.ones((1, 1, 1)))
    V = schur_transform(2, 5)
    assert all(len(lam) <= 2 for lam in V)


def _block_diag(blocks):
    size = sum(b.shape[0] for b in blocks)
    out = np.zeros((size, size))
    i = 0
    for b in blocks:
        out[i : i + b.shape[0], i : i + b.shape[0]] = b
        i += b.shape[0]
    return out
