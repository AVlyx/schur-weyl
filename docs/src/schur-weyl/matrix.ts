/**
 * Minimal dense-matrix helpers, standing in for the small slice of numpy
 * (`np.eye`, `np.zeros`, `@`, `np.linalg.det`) used by the Python library.
 * Not a general linear-algebra library — just enough for this port.
 */

export type Matrix = number[][];

export function zerosMatrix(rows: number, cols: number): Matrix {
  return Array.from({ length: rows }, () => new Array(cols).fill(0));
}

export function identityMatrix(n: number): Matrix {
  const m = zerosMatrix(n, n);
  for (let i = 0; i < n; i++) m[i][i] = 1;
  return m;
}

export function matMul(a: Matrix, b: Matrix): Matrix {
  const rows = a.length;
  const inner = b.length;
  const cols = inner === 0 ? 0 : b[0].length;
  const result = zerosMatrix(rows, cols);
  for (let i = 0; i < rows; i++) {
    for (let k = 0; k < inner; k++) {
      const aik = a[i][k];
      if (aik === 0) continue;
      for (let j = 0; j < cols; j++) {
        result[i][j] += aik * b[k][j];
      }
    }
  }
  return result;
}

/**
 * Determinant via Gaussian elimination with partial pivoting.
 *
 * Mirrors `np.linalg.det`: floating point, not exact for large integer
 * matrices. The Python library notes `schur_polynomial` is "not very
 * precise above len(lam) = 7" for the same reason.
 */
export function determinant(matrix: Matrix): number {
  const n = matrix.length;
  if (n === 0) return 1;
  const m = matrix.map((row) => [...row]);
  let det = 1;

  for (let col = 0; col < n; col++) {
    let pivotRow = col;
    let maxAbs = Math.abs(m[col][col]);
    for (let r = col + 1; r < n; r++) {
      const v = Math.abs(m[r][col]);
      if (v > maxAbs) {
        maxAbs = v;
        pivotRow = r;
      }
    }
    if (maxAbs === 0) return 0;
    if (pivotRow !== col) {
      [m[col], m[pivotRow]] = [m[pivotRow], m[col]];
      det = -det;
    }
    det *= m[col][col];
    for (let r = col + 1; r < n; r++) {
      const factor = m[r][col] / m[col][col];
      if (factor === 0) continue;
      for (let c = col; c < n; c++) {
        m[r][c] -= factor * m[col][c];
      }
    }
  }
  return det;
}
