/**
 * Ported from `src/schur_weyl/isotypic.py`.
 *
 * Note: the projector is a d^k x d^k dense matrix, so this is only
 * practical for small k and d (same limitation as the numpy original).
 */

import type { Partition } from './youngDiagrams';
import type { Permutation } from './symmetricGroup';
import { allPermutationsByCycleType, permutationInverse } from './symmetricGroup';
import { character } from './character';
import { dimSpecht } from './dimensions';
import { factorial } from './utils';

/**
 * Equivalent of:
 *   idx = np.arange(d**n).reshape((d,) * n)
 *   np.transpose(idx, axes=permutation_inverse(sigma)).reshape(-1)
 *
 * i.e. the permutation of flat indices of (C^d)^{\otimes n} induced by
 * permuting the tensor factors according to sigma.
 */
function permIndexMap(sigma: Permutation, d: number): number[] {
  const n = sigma.length;
  const axes = permutationInverse(sigma);
  const strides = new Array<number>(n);
  for (let m = 0; m < n; m++) strides[m] = d ** (n - 1 - m);
  const total = d ** n;

  const result = new Array<number>(total);
  for (let flat = 0; flat < total; flat++) {
    // decode flat -> digits, digit[0] most significant (row-major, like numpy)
    const digits = new Array<number>(n);
    let rem = flat;
    for (let m = 0; m < n; m++) {
      digits[m] = Math.floor(rem / strides[m]);
      rem %= strides[m];
    }
    // permuted[axes[k]] = digits[k]
    const permuted = new Array<number>(n);
    for (let k = 0; k < n; k++) permuted[axes[k]] = digits[k];
    let encoded = 0;
    for (let m = 0; m < n; m++) encoded += permuted[m] * strides[m];
    result[flat] = encoded;
  }
  return result;
}

/**
 * Return the isotypic projector acting on (C^d)^{\otimes k}, k = sum(lam).
 *
 * @param lam - the partition, k = sum(lam)
 * @param d - the basis dimension
 * @returns a dense d^k x d^k projector matrix
 */
export function isotypicProj(lam: Partition, d: number): number[][] {
  const k = lam.reduce((a, b) => a + b, 0);
  const dimension = d ** k;
  const proj: number[][] = Array.from({ length: dimension }, () => new Array<number>(dimension).fill(0));

  for (const { cycleType, perms } of allPermutationsByCycleType(k)) {
    const char = character(lam, cycleType);
    if (char === 0) continue;
    for (const perm of perms) {
      const mapped = permIndexMap(perm, d);
      for (let r = 0; r < dimension; r++) {
        proj[r][mapped[r]] += char;
      }
    }
  }

  const scale = dimSpecht(lam) / factorial(k);
  for (let r = 0; r < dimension; r++) {
    for (let c = 0; c < dimension; c++) {
      proj[r][c] *= scale;
    }
  }
  return proj;
}
