/**
 * What this application is willing to display, and why.
 *
 * Every fold below has a plausible energy term in `predict.ts`, and a plausible
 * energy term is exactly the problem: it produces a number, the number renders
 * as a bar in an ensemble, and a bar in an ensemble reads as a prediction. The
 * question a capability gate answers is not "can we compute something?" but
 * "is there evidence behind the thing we would draw?".
 *
 * The line is drawn where QuadCond's atlas ends. G-quadruplex and i-motif have
 * measured thermodynamics behind them -- 4,469 biophysical rows, seven
 * biophysically anchored heads. Duplex and hairpin have decades of published
 * nearest-neighbour parameters, which the local functions implement, so they are
 * shown and badged `heuristic`. AC-motif, cruciform and triplex have neither: no
 * head covers them and no row in the atlas constrains them, so what the app
 * would draw is an energy function's opinion with nothing checking it.
 *
 * RNA is a separate case and the strictest one. Every measurement in the atlas
 * is DNA. rG4s are effectively always parallel, more stable at the same ionic
 * strength, and the 2'-OH changes the loop energetics most of these features
 * encode. A DNA-trained head returns a confident number for an RNA sequence with
 * nothing behind it, and the service refuses rather than warns -- a warning
 * beside a number in a 3D viewer reads as a caveat on a real answer rather than
 * as the absence of one.
 */
import type { StructureKind } from "./types";

export interface Capability {
  enabled: boolean;
  /** Shown to the user wherever the fold would otherwise have appeared. */
  reason: string;
}

export const FOLD_CAPABILITIES: Record<StructureKind, Capability> = {
  "g-quadruplex": {
    enabled: true,
    reason: "Measured melting temperatures across 255 buffers; a condition-aware head.",
  },
  "i-motif": {
    enabled: true,
    reason:
      "Measured transitional pH and, for two constructs, a full condition response.",
  },
  duplex: {
    enabled: true,
    reason:
      "Published nearest-neighbour parameters. Shown as the competing state, badged heuristic.",
  },
  hairpin: {
    enabled: true,
    reason:
      "Published nearest-neighbour parameters. Shown as the competing state, badged heuristic.",
  },
  "ss-coil": {
    enabled: true,
    reason: "The unfolded reference state, at deltaG = 0 by definition.",
  },
  "ac-motif": {
    enabled: false,
    reason:
      "No QuadCond head covers the AC-motif and no row in the atlas constrains it. " +
      "The local energy term would produce a number with nothing behind it.",
  },
  cruciform: {
    enabled: false,
    reason:
      "No QuadCond head covers cruciform extrusion. Predicting it needs superhelical " +
      "density, which this tool does not model.",
  },
  triplex: {
    enabled: false,
    reason:
      "No QuadCond head covers triplex formation, and it requires a third strand or a " +
      "mirror-repeat context the single-sequence input does not supply.",
  },
};

/** Protein binders were stored as truncated sequences; the panel is off until fixed. */
export const PROTEIN_BINDERS: Capability = {
  enabled: false,
  reason:
    "The bundled binder records are unreliable -- BLM (UniProt P54132) was stored as " +
    "249 residues against a true length of 1,417 -- and binding was scored by motif " +
    "class rather than by any measurement. Nothing here is evidence about binding.",
};

export const RNA_SUPPORT: Capability = {
  enabled: false,
  reason:
    "Every measurement in the QuadCond atlas is DNA. rG4s are effectively always " +
    "parallel, more stable at the same ionic strength, and the 2'-OH changes the loop " +
    "energetics these features encode. The service refuses RNA rather than answering it.",
};

export function foldEnabled(kind: StructureKind): boolean {
  return FOLD_CAPABILITIES[kind]?.enabled ?? false;
}

/** The folds withheld, for the interface to list rather than silently omit. */
export function disabledFolds(): Array<{ kind: StructureKind; reason: string }> {
  return (Object.keys(FOLD_CAPABILITIES) as StructureKind[])
    .filter((k) => !FOLD_CAPABILITIES[k].enabled)
    .map((k) => ({ kind: k, reason: FOLD_CAPABILITIES[k].reason }));
}
