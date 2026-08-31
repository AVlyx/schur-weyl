/**
 * Ported from `src/schur_weyl/tableaux.py`.
 */

import type { Partition } from './youngDiagrams';
import { removableCorners } from './youngDiagrams';

/** A (semi-)standard Young tableau: rows of entries, ragged like the shape. */
export type Tableau = number[][];

function cloneTableau(t: Tableau): Tableau {
  return t.map((row) => [...row]);
}

/**
 * Generate all standard Young tableaux of a given shape.
 *
 * @example
 * [...standardYoungTableaux([])] // [[]]
 * [...standardYoungTableaux([2, 1])] // [[[1,2],[3]], [[1,3],[2]]]
 * [...standardYoungTableaux([3, 2, 1])].length // 16
 */
export function* standardYoungTableaux(lam: Partition): Generator<Tableau, void, void> {
  if (lam.length === 0) {
    yield [];
    return;
  }

  const ret: number[][] = lam.map((li) => new Array<number>(li).fill(-1));
  ret[0][0] = 1;

  function* sytRec(smallLam: number[], d: number): Generator<Tableau, void, void> {
    if (smallLam.length === 1 && smallLam[0] === 1) {
      yield cloneTableau(ret);
      return;
    }
    for (const [i, j] of removableCorners(smallLam)) {
      if (j === 0) {
        smallLam.splice(i, 1);
        ret[i][j] = d;
        yield* sytRec(smallLam, d - 1);
        smallLam.splice(i, 0, 1);
      } else {
        smallLam[i] -= 1;
        ret[i][j] = d;
        yield* sytRec(smallLam, d - 1);
        smallLam[i] += 1;
      }
    }
  }

  yield* sytRec([...lam], lam.reduce((a, b) => a + b, 0));
}

/**
 * Generate all semi-standard Young tableaux of a given shape and alphabet.
 *
 * @param d - the alphabet size
 *
 * @example
 * [...semiStandardYoungTableau([2, 1], 3)]
 * // [[[1,1],[2]], [[1,2],[2]], [[1,3],[2]], [[1,1],[3]], [[1,2],[3]],
 * //  [[1,3],[3]], [[2,2],[3]], [[2,3],[3]]]
 */
export function* semiStandardYoungTableau(lam: Partition, d: number): Generator<Tableau, void, void> {
  if (lam.length === 0) {
    yield [];
    return;
  }
  const ret: number[][] = lam.map((li) => new Array<number>(li).fill(-1));
  const maxv: number[][] = lam.map((li) => new Array<number>(li).fill(d + 1));

  function* fill(i: number, j: number, minv: number): Generator<Tableau, void, void> {
    if (j >= lam[i]) {
      if (i === 0) {
        yield cloneTableau(ret);
        return;
      }
      yield* fill(i - 1, 0, i);
      return;
    }
    for (let cellVal = minv; cellVal < maxv[i][j]; cellVal++) {
      ret[i][j] = cellVal;
      if (i !== 0) maxv[i - 1][j] = cellVal;
      yield* fill(i, j + 1, cellVal);
    }
  }

  yield* fill(lam.length - 1, 0, lam.length);
}

/**
 * Generate all semi-standard Young tableaux with given content.
 *
 * @param content - should satisfy sum(lam) === sum(content)
 *
 * @example
 * [...semiStandardYoungTableauByContent([2, 1], [2, 1])] // [[[1,1],[2]]]
 * [...semiStandardYoungTableauByContent([2, 1], [1, 1, 1])]
 * // [[[1,3],[2]], [[1,2],[3]]]
 * [...semiStandardYoungTableauByContent([], [])] // [[]]
 */
export function* semiStandardYoungTableauByContent(
  lam: Partition,
  content: number[],
): Generator<Tableau, void, void> {
  const lamSum = lam.reduce((a, b) => a + b, 0);
  const contentSum = content.reduce((a, b) => a + b, 0);
  if (lamSum !== contentSum) return;
  if (lam.length === 0) {
    yield [];
    return;
  }

  const d = content.length;
  const contentCp = [...content];
  const ret: number[][] = lam.map((li) => new Array<number>(li).fill(-1));
  const maxv: number[][] = lam.map((li) => new Array<number>(li).fill(d + 1));

  function* fill(i: number, j: number, minv: number): Generator<Tableau, void, void> {
    if (j >= lam[i]) {
      if (i === 0) {
        yield cloneTableau(ret);
        return;
      }
      yield* fill(i - 1, 0, i);
      return;
    }
    for (let cellVal = minv; cellVal < maxv[i][j]; cellVal++) {
      if (contentCp[cellVal - 1] <= 0) continue;
      ret[i][j] = cellVal;
      if (i !== 0) maxv[i - 1][j] = cellVal;
      contentCp[cellVal - 1] -= 1;
      yield* fill(i, j + 1, cellVal);
      contentCp[cellVal - 1] += 1;
    }
  }

  yield* fill(lam.length - 1, 0, lam.length);
}

// TODO implement with a better algorithm (horizontal strips)
/**
 * Compute the Kostka number K_{lam, mu}.
 *
 * @example
 * kostka([2, 1], [1, 1, 1]) // 2
 */
export function kostka(lam: Partition, mu: Partition): number {
  return [...semiStandardYoungTableauByContent(lam, mu)].length;
}

/**
 * The reverse reading word of a tableau: right-to-left within each row,
 * bottom row to top row.
 *
 * @example
 * reverseReadingWord([[1,1,3],[2,3],[4]]) // [4,3,2,3,1,1]
 * reverseReadingWord([[1,2],[3]]) // [3,2,1]
 * reverseReadingWord([]) // []
 */
export function reverseReadingWord(tableau: Tableau): number[] {
  const result: number[] = [];
  for (let r = tableau.length - 1; r >= 0; r--) {
    const row = tableau[r];
    for (let c = row.length - 1; c >= 0; c--) result.push(row[c]);
  }
  return result;
}
