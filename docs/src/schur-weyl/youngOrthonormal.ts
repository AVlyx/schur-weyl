/**
 * Ported from `src/schur_weyl/young_orthonormal.py`.
 * Young's orthogonal form: explicit matrices for the irreps of S_n.
 *
 * Conventions (consistent with the rest of the project):
 *   - cells are 0-indexed, content(i, j) = j - i
 *   - tableau ENTRIES are 1-indexed, i.e. 1..n
 *   - permutations are 0-indexed one-line arrays, compose(s, t)(i) = s(t(i))
 *   - the adjacent transposition s_k swaps the ENTRIES k and k+1, for
 *     1 <= k <= n-1. It equals the 0-indexed permutation swapping k-1 and k.
 */

import type { Partition } from './youngDiagrams';
import type { Tableau } from './tableaux';
import type { Permutation } from './symmetricGroup';
import type { Matrix } from './matrix';
import { standardYoungTableaux } from './tableaux';
import { dimSpecht } from './dimensions';
import { identityMatrix, matMul } from './matrix';

function stripZeros(lam: Partition): Partition {
  return lam.filter((x) => x > 0);
}

// ------------------------------------------------------------------ tableaux

const sytBasisCache = new Map<string, Tableau[]>();

export function sytBasis(lam: Partition): Tableau[] {
  const clean = stripZeros(lam);
  const key = JSON.stringify(clean);
  let basis = sytBasisCache.get(key);
  if (!basis) {
    basis = [...standardYoungTableaux(clean)];
    sytBasisCache.set(key, basis);
  }
  return basis;
}

/** {entry: [i, j]} for a tableau. */
function positions(T: Tableau): Map<number, [number, number]> {
  const pos = new Map<number, [number, number]>();
  T.forEach((row, i) => {
    row.forEach((v, j) => pos.set(v, [i, j]));
  });
  return pos;
}

/** The tableau with entries k and k+1 exchanged (may not be standard). */
function swapEntries(T: Tableau, k: number): Tableau {
  return T.map((row) => row.map((v) => (v === k ? k + 1 : v === k + 1 ? k : v)));
}

function tableauKey(T: Tableau): string {
  return JSON.stringify(T);
}

// ---------------------------------------------------------- adjacent generator

/** cont(box holding k+1) - cont(box holding k), with content = col - row. */
export function axialDistance(T: Tableau, k: number): number {
  const pos = positions(T);
  const p1 = pos.get(k);
  const p2 = pos.get(k + 1);
  if (!p1 || !p2) throw new Error(`entries ${k}, ${k + 1} not found in tableau`);
  const [i1, j1] = p1;
  const [i2, j2] = p2;
  return j2 - i2 - (j1 - i1);
}

const generatorCache = new Map<string, Matrix>();

/**
 * The matrix of rho^lam(s_k), where s_k swaps the entries k and k+1.
 *
 * On the span of {T, s_k T}, with r = 1 / axialDistance(T, k):
 *   [  r              sqrt(1 - r^2) ]
 *   [  sqrt(1 - r^2)  -r            ]
 *
 * If k and k+1 share a row the distance is +1 (so r = 1, eigenvalue +1);
 * if they share a column it is -1 (eigenvalue -1). In both cases s_k T is
 * not standard and the off-diagonal term vanishes on its own.
 */
export function youngOrthogonalGenerator(lam: Partition, k: number): Matrix {
  const clean = stripZeros(lam);
  const key = JSON.stringify([clean, k]);
  const cached = generatorCache.get(key);
  if (cached) return cached;

  const basis = sytBasis(clean);
  const index = new Map<string, number>();
  basis.forEach((T, a) => index.set(tableauKey(T), a));
  const f = basis.length;
  const M = Array.from({ length: f }, () => new Array<number>(f).fill(0));

  basis.forEach((T, b) => {
    const delta = axialDistance(T, k);
    const r = 1 / delta;
    M[b][b] = r;
    if (Math.abs(delta) === 1) {
      // same row (+1) or same column (-1): no mixing partner
      return;
    }
    const a = index.get(tableauKey(swapEntries(T, k)));
    if (a === undefined) throw new Error('expected swapped tableau in basis');
    M[a][b] = Math.sqrt(1 - r * r);
  });

  generatorCache.set(key, M);
  return M;
}

// ------------------------------------------------------------ general element

/** [i_1, ..., i_m] (0-indexed) with sigma == s_{i_1} s_{i_2} ... s_{i_m}. */
function adjacentFactorization(sigma: Permutation): number[] {
  const arr = [...sigma];
  const n = arr.length;
  const right: number[] = [];
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n - 1 - i; j++) {
      if (arr[j] > arr[j + 1]) {
        [arr[j], arr[j + 1]] = [arr[j + 1], arr[j]];
        right.push(j);
      }
    }
  }
  return right.reverse();
}

/**
 * rho^lam(sigma) for a 0-indexed one-line permutation sigma.
 * Factors sigma into adjacent transpositions and multiplies the generators.
 */
export function youngOrthogonal(lam: Partition, sigma: Permutation): Matrix {
  const clean = stripZeros(lam);
  const f = dimSpecht(clean);
  let M = identityMatrix(f);
  for (const i of adjacentFactorization(sigma)) {
    M = matMul(M, youngOrthogonalGenerator(clean, i + 1)); // 0-indexed i -> entry i+1
  }
  return M;
}

// --------------------------------------------------------- Jucys-Murphy

/**
 * rho^lam(X_k) with X_k = sum_{i<k} (i k), 1-indexed: X_k = sum over the
 * transpositions swapping entries i and k for i = 1..k-1.
 *
 * In Young's orthogonal basis this is DIAGONAL, with the entry for tableau T
 * equal to the content of the box holding k.
 */
export function jucysMurphy(lam: Partition, k: number): Matrix {
  const clean = stripZeros(lam);
  const n = clean.reduce((a, b) => a + b, 0);
  const f = dimSpecht(clean);
  const total = Array.from({ length: f }, () => new Array<number>(f).fill(0));
  for (let i = 1; i < k; i++) {
    // 0-indexed transposition swapping i-1 and k-1
    const p = Array.from({ length: n }, (_, idx) => idx);
    [p[i - 1], p[k - 1]] = [p[k - 1], p[i - 1]];
    const Y = youngOrthogonal(clean, p);
    for (let r = 0; r < f; r++) {
      for (let c = 0; c < f; c++) total[r][c] += Y[r][c];
    }
  }
  return total;
}

/**
 * The predicted diagonal of jucysMurphy(lam, k): content of the box holding
 * the entry k, one value per standard tableau.
 */
export function jucysMurphyContents(lam: Partition, k: number): number[] {
  const clean = stripZeros(lam);
  const out: number[] = [];
  for (const T of sytBasis(clean)) {
    const pos = positions(T).get(k);
    if (!pos) throw new Error(`entry ${k} not found in tableau`);
    const [i, j] = pos;
    out.push(j - i);
  }
  return out;
}
