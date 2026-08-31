/** Small shared helpers used across the schur-weyl port. */

/**
 * n! for non-negative integers.
 *
 * JS numbers are IEEE-754 doubles, so this is exact only up to n = 18
 * (18! < 2^53). The Python original uses arbitrary-precision integers;
 * this port trades that off for a dependency-free implementation, which
 * is fine for the partition sizes this library is used with in the docs.
 */
export function factorial(n: number): number {
  let result = 1;
  for (let i = 2; i <= n; i++) result *= i;
  return result;
}
