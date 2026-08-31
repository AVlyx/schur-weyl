/**
 * TypeScript port of the `schur_weyl` Python library
 * (https://github.com/AVlyx/schur-weyl, `src/schur_weyl/`).
 *
 * The re-exports below mirror `src/schur_weyl/__init__.py`; the individual
 * modules also export a few extra helpers (permutations, characters, Young's
 * orthogonal form, the Schur-Weyl measure) that aren't part of the Python
 * package's public surface but are ported here for the docs site to use.
 */

export type { Partition } from './youngDiagrams';
export {
  partitions,
  partitionConjugate,
  hookLength,
  cells,
  majorizes,
  addableCorners,
  removableCorners,
  isInPartition,
} from './youngDiagrams';

export { dimSpecht, dimWeyl } from './dimensions';

export type { Tableau } from './tableaux';
export {
  standardYoungTableaux,
  semiStandardYoungTableau,
  semiStandardYoungTableauByContent,
  kostka,
  reverseReadingWord,
} from './tableaux';

export { schurPolynomial, powerSum, h0ToK } from './symmetricFunctions';

export { isotypicProj } from './isotypic';

// --- extras beyond src/schur_weyl/__init__.py's public surface ---

export type { Permutation, CycleClass } from './symmetricGroup';
export {
  permutationIdentity,
  permutationCompose,
  permutationInverse,
  allPermutations,
  allPermutationsByCycleType,
  permutationCycleType,
  permutationConjugacyClassSize,
} from './symmetricGroup';

export type { CharacterTable } from './character';
export { character, characterTable, borderStrips } from './character';

export type { SchurWeylEntry } from './swMeasure';
export { schurWeylMeasure } from './swMeasure';

export {
  sytBasis,
  axialDistance,
  youngOrthogonalGenerator,
  youngOrthogonal,
  jucysMurphy,
  jucysMurphyContents,
} from './youngOrthonormal';

export type { Matrix } from './matrix';
