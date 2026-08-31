/**
 * Ported from `src/schur_weyl/symmetric_group.py`.
 *
 * A permutation of {0, ..., n-1} is represented as its 0-indexed one-line
 * array, e.g. `[1, 2, 0]` sends 0 -> 1, 1 -> 2, 2 -> 0.
 */

import { factorial } from './utils';

export type Permutation = number[];

/**
 * The identity permutation on {0, ..., n-1}.
 *
 * @example
 * permutationIdentity(4) // [0, 1, 2, 3]
 */
export function permutationIdentity(n: number): Permutation {
  return Array.from({ length: n }, (_, i) => i);
}

/**
 * Multiply two permutations: (pi * sigma)(i) = pi(sigma(i)).
 *
 * @example
 * permutationCompose([1, 0, 2], [0, 2, 1]) // [1, 2, 0]
 * permutationCompose([0, 1, 2], [2, 0, 1]) // [2, 0, 1]
 */
export function permutationCompose(pi: Permutation, sigma: Permutation): Permutation {
  return sigma.map((_, i) => pi[sigma[i]]);
}

/**
 * The inverse permutation.
 *
 * @example
 * permutationInverse([1, 2, 0]) // [2, 0, 1]
 * permutationInverse([0, 1, 2]) // [0, 1, 2]
 */
export function permutationInverse(sigma: Permutation): Permutation {
  const result = new Array<number>(sigma.length).fill(0);
  sigma.forEach((s, i) => {
    result[s] = i;
  });
  return result;
}

/** All permutations of {0, ..., n-1}. */
export function* allPermutations(n: number): Generator<Permutation, void, void> {
  yield* permuteRange(Array.from({ length: n }, (_, i) => i));
}

function* permuteRange(arr: number[]): Generator<Permutation, void, void> {
  if (arr.length <= 1) {
    yield [...arr];
    return;
  }
  for (let i = 0; i < arr.length; i++) {
    const rest = [...arr.slice(0, i), ...arr.slice(i + 1)];
    for (const perm of permuteRange(rest)) {
      yield [arr[i], ...perm];
    }
  }
}

/**
 * Give the cycle type of a permutation: the sizes of its disjoint cycles,
 * sorted in decreasing order.
 *
 * @example
 * permutationCycleType([1, 2, 0, 4, 3]) // [3, 2]
 * permutationCycleType([0, 1, 2]) // [1, 1, 1]
 * permutationCycleType([1, 0, 2]) // [2, 1]
 */
export function permutationCycleType(perm: Permutation): number[] {
  const seen = new Array<boolean>(perm.length).fill(false);
  const cycleLens: number[] = [];
  for (let i = 0; i < perm.length; i++) {
    if (seen[i]) continue;
    seen[i] = true;
    let count = 1;
    let j = perm[i];
    while (!seen[j]) {
      seen[j] = true;
      j = perm[j];
      count += 1;
    }
    cycleLens.push(count);
  }
  return cycleLens.sort((a, b) => b - a);
}

/** All permutations of {0, ..., n-1}, grouped by cycle type. */
export interface CycleClass {
  cycleType: number[];
  perms: Permutation[];
}

export function allPermutationsByCycleType(n: number): CycleClass[] {
  const byKey = new Map<string, CycleClass>();
  for (const p of allPermutations(n)) {
    const cycleType = permutationCycleType(p);
    const key = cycleType.join(',');
    let entry = byKey.get(key);
    if (!entry) {
      entry = { cycleType, perms: [] };
      byKey.set(key, entry);
    }
    entry.perms.push(p);
  }
  return [...byKey.values()];
}

/**
 * The number of permutations of the same cycle type as cycleType (the size
 * of its conjugacy class in S_n, n = sum(cycleType)).
 *
 * @example
 * permutationConjugacyClassSize([1, 1, 1, 1]) // 1
 * permutationConjugacyClassSize([2, 1, 1]) // 6
 */
export function permutationConjugacyClassSize(cycleType: number[]): number {
  const counts = new Map<number, number>();
  for (const c of cycleType) counts.set(c, (counts.get(c) ?? 0) + 1);
  let z = 1;
  for (const [cycleLen, cycleCount] of counts) {
    z *= cycleLen ** cycleCount * factorial(cycleCount);
  }
  const n = cycleType.reduce((a, b) => a + b, 0);
  return Math.round(factorial(n) / z);
}
