/**
 * Ported from `src/schur_weyl/character.py`.
 * Characters of the symmetric group via the Murnaghan-Nakayama rule.
 */

import type { Partition } from './youngDiagrams';
import { partitions } from './youngDiagrams';

// ------------------------------------------------------------- border strips

/**
 * Yield [newShape, height] for every border strip of size k removable from
 * lam.
 *
 * Uses beta numbers (first-column hook lengths). With L = lam.length, set
 *   beta_i = lam_i + (L - 1 - i),
 * a strictly decreasing sequence of L distinct non-negative integers.
 * Removing a border strip of size k is exactly: pick i with
 * beta_i - k >= 0 and beta_i - k not already a beta value, then replace
 * beta_i by beta_i - k. The height of the strip is the number of beta
 * values jumped over.
 */
export function* borderStrips(lam: Partition, k: number): Generator<[Partition, number], void, void> {
  const L = lam.length;
  if (L === 0) return;
  const beta = lam.map((v, i) => v + (L - 1 - i));
  const bset = new Set(beta);
  for (let i = 0; i < L; i++) {
    const target = beta[i] - k;
    if (target < 0 || bset.has(target)) continue;
    let height = 0;
    for (const b of beta) {
      if (target < b && b < beta[i]) height++;
    }
    const newBeta = [...beta.slice(0, i), target, ...beta.slice(i + 1)].sort((a, b) => b - a);
    const newLam = newBeta.map((b, j) => b - (L - 1 - j));
    while (newLam.length > 0 && newLam[newLam.length - 1] === 0) newLam.pop();
    yield [newLam, height];
  }
}

// ------------------------------------------------------- Murnaghan-Nakayama

const characterCache = new Map<string, number>();

function characterMemo(lam: Partition, mu: Partition): number {
  const key = JSON.stringify([lam, mu]);
  const cached = characterCache.get(key);
  if (cached !== undefined) return cached;

  let result: number;
  if (mu.length === 0) {
    result = lam.length === 0 ? 1 : 0;
  } else {
    const k = mu[0];
    const rest = mu.slice(1);
    let total = 0;
    for (const [newLam, height] of borderStrips(lam, k)) {
      total += (height % 2 === 0 ? 1 : -1) * characterMemo(newLam, rest);
    }
    result = total;
  }

  characterCache.set(key, result);
  return result;
}

/**
 * chi^lam(cycleType): the irreducible character of S_n indexed by the
 * partition lam, evaluated on the conjugacy class of cycle type cycleType.
 *
 * @example
 * character([3], [1, 1, 1]) // trivial rep -> 1
 * character([1, 1, 1], [2, 1]) // sign rep at a transposition -> -1
 * character([2, 1], [1, 1, 1]) // dimension of the standard rep -> 2
 * character([2, 1], [3]) // -1
 * character([2, 1], [2, 1]) // 0
 */
export function character(lam: Partition, cycleType: Partition): number {
  const ct = cycleType
    .filter((x) => x > 0)
    .sort((a, b) => b - a);
  const lamSum = lam.reduce((a, b) => a + b, 0);
  const ctSum = ct.reduce((a, b) => a + b, 0);
  if (lamSum !== ctSum) {
    throw new Error(`size mismatch: |${JSON.stringify(lam)}| != |${JSON.stringify(ct)}|`);
  }
  return characterMemo(lam, ct);
}

export interface CharacterTable {
  rows: Partition[];
  cols: Partition[];
  table: number[][];
}

/**
 * (rows, cols, table) with rows/cols the partitions of n in the order
 * produced by `partitions(n)`, and table[a][b] = chi^{rows[a]}(cols[b]).
 */
export function characterTable(n: number): CharacterTable {
  const ps = [...partitions(n)];
  const table = ps.map((lam) => ps.map((mu) => character(lam, mu)));
  return { rows: ps, cols: ps, table };
}
