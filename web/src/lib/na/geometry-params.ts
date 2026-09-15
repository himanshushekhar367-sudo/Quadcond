/**
 * Helical parameters for the folds this app draws.
 *
 * These are separated from the builders because the previous version had no
 * such file, and that is exactly how the bug happened: `geometry.ts` carried
 * one set of constants under the comment "Display-scale coordinates" and used
 * them for every fold. The values it used — 3.38 Å rise, ~34° twist, ~9 Å
 * radius — are correct B-DNA numbers. They are wrong for a G-quadruplex, which
 * is four strands around stacked tetrads, and wrong by a factor of two for an
 * i-motif, whose two duplexes intercalate.
 *
 * Nothing here is fitted to a particular PDB entry. They are round numbers from
 * the consensus ranges, which is the honest precision for an idealised model.
 * Every structure built from them is labelled a MODEL, and the PDB export says
 * so in its REMARK records.
 *
 * The Python reference implementation in `quadcond/structure/geometry.py`
 * measures what it builds — tetrad planarity, inter-tetrad rise, four-fold
 * occupancy, backbone step lengths and a hard-sphere clash test — and these
 * constants are the ones that pass it.
 */

/** B-form duplex. Unchanged: the duplex builder was already right. */
export const DUPLEX = {
  rise: 3.38,
  twistDeg: 34.3,
  c1Radius: 9.0,
} as const;

/**
 * G-quadruplex. Four columns at 90°, stacked tetrads.
 *
 * `riseAngstrom` is the stacking distance between successive G-tetrads;
 * `twistDeg` is the right-handed rotation between them. `c1Radius` places the
 * sugar C1′ atoms relative to the four-fold axis and `baseRadius` the guanine
 * ring centres, whose O6 atoms coordinate the axial cation — which is why the
 * base ring sits far inside the backbone.
 */
export const G4 = {
  rise: 3.3,
  twistDeg: 30.0,
  c1Radius: 9.7,
  baseRadius: 4.6,
} as const;

/**
 * i-motif. Two parallel-stranded duplexes intercalated with antiparallel
 * polarity.
 *
 * The two rises are the point. Successive C:C⁺ pairs *of one duplex* are
 * `duplexRise` apart because the partner duplex's pairs sit between them; the
 * intercalated stack as a whole advances `pairRise` per step. A model built on
 * a single duplex rise gets the height of the fold wrong by about two-fold and
 * has nowhere to put the intercalated pairs.
 */
export const IMOTIF = {
  pairRise: 3.1,
  duplexRise: 6.2,
  twistDeg: 12.0,
  c1Radius: 5.4,
  pairSeparation: 8.6,
} as const;

/** Typical P–P separation along a strand, used to route loops and flanks. */
export const BACKBONE_STEP = 6.0;

/** Hard-sphere floor: no two placed C1′ atoms may come closer than this. */
export const MIN_CONTACT = 2.8;

export const MODEL_DISCLAIMER =
  "Idealised model geometry at published helical parameters. Not a refined " +
  "structure, not a prediction of atomic coordinates, and not fitted to any " +
  "deposited entry.";
