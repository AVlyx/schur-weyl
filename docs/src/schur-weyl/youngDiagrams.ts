/**
 * Ported from `src/schur_weyl/young_diagrams.py`.
 *
 * A partition is represented as a weakly-decreasing array of positive
 * integers (no trailing zeros), e.g. `[4, 3, 1]` for the Python tuple
 * `(4, 3, 1)`.
 */

export type Partition = number[];

/**
 * Generate the partitions of k: all non-increasing arrays of positive
 * integers summing to k.
 *
 * @param k - the size
 * @param maxPart - the maximum size of any part. Defaults to k.
 * @param maxHeight - the maximum number of parts. Defaults to k.
 *
 * @example
 * [...partitions(4)] // [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
 * [...partitions(5, 2)] // [[2, 2, 1], [2, 1, 1, 1], [1, 1, 1, 1, 1]]
 * [...partitions(5, 2, 4)] // [[2, 2, 1], [2, 1, 1, 1]]
 */
export function* partitions(
  k: number,
  maxPart?: number,
  maxHeight?: number,
): Generator<Partition, void, void> {
  const mp = maxPart ?? k;
  const mh = maxHeight ?? k;
  if (k === 0) {
    yield [];
    return;
  }
  if (mh === 0) {
    return;
  }

  for (let first = Math.min(k, mp); first >= 1; first--) {
    for (const rest of partitions(k - first, first, mh - 1)) {
      yield [first, ...rest];
    }
  }
}

/**
 * Take the conjugate of a partition: flip/transpose a Young diagram across
 * its main diagonal.
 *
 * @example
 * partitionConjugate([4, 3, 1]) // [3, 2, 2, 1]
 * partitionConjugate([]) // []
 */
export function partitionConjugate(lam: Partition): Partition {
  if (lam.length === 0) return [];
  const result = new Array<number>(lam[0]).fill(0);
  for (const p of lam) {
    for (let j = 0; j < p; j++) result[j] += 1;
  }
  return result;
}

/**
 * The hook length at cell (i, j): (arm + 1) + leg.
 *
 * @param i - the row index (0-indexed)
 * @param j - the column index (0-indexed)
 *
 * @example
 * hookLength([4, 3, 1], 0, 0) // 6
 * hookLength([4, 3, 3, 2], 2, 1) // 3
 */
export function hookLength(lam: Partition, i: number, j: number): number {
  let leg = 0;
  for (let l = i + 1; l < lam.length; l++) {
    if (lam[l] > j) leg++;
  }
  return lam[i] - j + leg;
}

/**
 * All the cell indices of a partition (0-indexed).
 *
 * @example
 * [...cells([3, 2, 1])] // [[0,0],[0,1],[0,2],[1,0],[1,1],[2,0]]
 */
export function* cells(lam: Partition): Generator<[number, number], void, void> {
  for (let i = 0; i < lam.length; i++) {
    for (let j = 0; j < lam[i]; j++) {
      yield [i, j];
    }
  }
}

/**
 * Check if lam majorizes nu. Assumes lam and nu have the same size,
 * otherwise the result is undefined.
 *
 * @example
 * majorizes([3, 3], [4, 1, 1]) // false
 * majorizes([4, 1, 1], [3, 3]) // false
 * majorizes([3, 1, 1, 1], [2, 2, 1, 1]) // true
 * majorizes([2, 1], [2, 1]) // true
 */
export function majorizes(lam: Partition, nu: Partition): boolean {
  if (nu.length < lam.length) return false;
  let cLam = 0;
  let cNu = 0;
  for (let i = 0; i < lam.length; i++) {
    cLam += lam[i];
    cNu += nu[i];
    if (cLam < cNu) return false;
  }
  return true;
}

/**
 * The cells that can be added to lam to still obtain a valid partition.
 *
 * @example
 * [...addableCorners([2, 2, 1])] // [[0,2],[2,1],[3,0]]
 */
export function* addableCorners(lam: Partition): Generator<[number, number], void, void> {
  if (lam.length === 0) {
    yield [0, 0];
    return;
  }
  yield [0, lam[0]];
  for (let i = 1; i < lam.length; i++) {
    if (lam[i] !== lam[i - 1]) yield [i, lam[i]];
  }
  yield [lam.length, 0];
}

/**
 * The cells that can be removed from lam to still obtain a valid partition.
 *
 * @example
 * [...removableCorners([2, 2, 1])] // [[2,0],[1,1]]
 */
export function* removableCorners(lam: Partition): Generator<[number, number], void, void> {
  if (lam.length === 0) return;
  yield [lam.length - 1, lam[lam.length - 1] - 1];
  for (let i = lam.length - 2; i >= 0; i--) {
    if (lam[i] !== lam[i + 1]) yield [i, lam[i] - 1];
  }
}

/**
 * True if the cell (i, j) is in the partition lam.
 *
 * @example
 * isInPartition([2, 1], 1, 1) // false
 * isInPartition([2, 1], 1, 0) // true
 */
export function isInPartition(lam: Partition, i: number, j: number): boolean {
  return lam.length > i && lam[i] > j;
}
