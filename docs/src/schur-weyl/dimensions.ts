/**
 * Ported from `src/schur_weyl/dimensions.py`.
 */

import type { Partition } from './youngDiagrams';
import { cells, hookLength } from './youngDiagrams';
import { factorial } from './utils';

/**
 * The dimension of the Specht module S^lam: the number of SYT of shape lam.
 * f^lam = n! / prod(hooks).
 *
 * @example
 * dimSpecht([2, 1]) // 2
 * dimSpecht([3, 2, 1]) // 16
 */
export function dimSpecht(lam: Partition): number {
  const n = lam.reduce((a, b) => a + b, 0);
  let hooks = 1;
  for (const [i, j] of cells(lam)) {
    hooks *= hookLength(lam, i, j);
  }
  return Math.round(factorial(n) / hooks);
}

/**
 * The dimension of the Weyl module V_lam for GL_d: the number of SSYT with
 * entries from 1 to d.
 *
 * Returns 0 if lam has more than d rows (a cell then has d + content = 0).
 *
 * @param d - the alphabet size
 *
 * @example
 * dimWeyl([2, 1], 2) // 2
 * dimWeyl([2, 1], 3) // 8
 * dimWeyl([2], 3) // 6
 */
export function dimWeyl(lam: Partition, d: number): number {
  let numerator = 1;
  let denominator = 1;
  for (const [i, j] of cells(lam)) {
    numerator *= d + j - i;
    denominator *= hookLength(lam, i, j);
  }
  return Math.round(numerator / denominator);
}
