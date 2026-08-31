/**
 * Ported from `src/schur_weyl/sw_measure.py`.
 */

import type { Partition } from './youngDiagrams';
import { partitions } from './youngDiagrams';
import { dimSpecht } from './dimensions';
import { schurPolynomial } from './symmetricFunctions';

export interface SchurWeylEntry {
  partition: Partition;
  probability: number;
}

/**
 * The probabilities tr(Pi^lam rho^{\otimes k}) where rho has the given
 * spectrum, one entry per partition of k.
 *
 * @param spectrum - the spectrum of rho
 * @param k - the number of copies of rho
 */
export function schurWeylMeasure(spectrum: number[], k: number): SchurWeylEntry[] {
  const ret: SchurWeylEntry[] = [];
  for (const lam of partitions(k)) {
    ret.push({ partition: lam, probability: dimSpecht(lam) * schurPolynomial(lam, spectrum) });
  }
  return ret;
}
