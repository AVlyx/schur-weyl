/**
 * Ported from `src/schur_weyl/symmetric_functions.py`.
 */

import type { Partition } from './youngDiagrams';
import { determinant } from './matrix';

/**
 * The Schur polynomial of a partition evaluated on input xs.
 * Not very precise above lam.length = 7 (see {@link determinant}).
 */
export function schurPolynomial(lam: Partition, xs: number[]): number {
  return jacobiTrudi(lam, xs);
}

/** Not very precise above lam.length = 7. Could use a different det formula. */
function jacobiTrudi(lam: Partition, xs: number[]): number {
  if (lam.length === 0) return 1.0;
  const ell = lam.length;
  const hk = h0ToK(ell + lam[0] - 1, xs);
  const mat: number[][] = [];
  for (let i = 0; i < ell; i++) {
    const row: number[] = [];
    for (let j = 0; j < ell; j++) {
      const idx = lam[i] - i + j;
      row.push(idx >= 0 ? hk[idx] : 0.0);
    }
    mat.push(row);
  }
  return determinant(mat);
}

/**
 * Completely homogeneous symmetric polynomials h_0 through h_k.
 *
 * @returns an array of length k + 1, [h_0, ..., h_k]
 */
export function h0ToK(k: number, xs: number[]): number[] {
  const h: number[] = [1.0, ...new Array<number>(k).fill(0.0)];
  for (const xi of xs) {
    for (let i = 1; i <= k; i++) {
      h[i] += xi * h[i - 1];
    }
  }
  return h;
}

/**
 * p_mu(x) = prod_i (sum_j x_j^{mu_i}).
 *
 * @example
 * powerSum([1], [2, 3]) // p_1 = 2 + 3 = 5
 * powerSum([2], [2, 3]) // p_2 = 4 + 9 = 13
 * powerSum([2, 1], [2, 3]) // p_2 * p_1 = 13 * 5 = 65
 * powerSum([], [2, 3]) // 1
 */
export function powerSum(mu: Partition, xs: number[]): number {
  let result = 1;
  for (const part of mu) {
    result *= xs.reduce((acc, xj) => acc + xj ** part, 0);
  }
  return result;
}
