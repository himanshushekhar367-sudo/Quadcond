export type Polymer = "DNA" | "RNA";

export type StructureKind =
  | "duplex"
  | "hairpin"
  | "g-quadruplex"
  | "i-motif"
  | "ac-motif"
  | "cruciform"
  | "triplex"
  | "ss-coil";

/**
 * The buffer, in the eight fields QuadCond's atlas actually records.
 *
 * The previous shape had a single `monovalent` slider that the client mapped to
 * K+. That is not a display simplification -- potassium and sodium have opposite
 * effects on i-motif stability and a ~13 degC difference on G4 melting at the
 * same concentration, which the `g4_tm` head separates out-of-fold at 94% sign
 * agreement. Collapsing them threw away the single best-established fact in G4
 * buffer chemistry and silently answered a question the user did not ask.
 *
 * Field names and units match `quadcond.conditions.Condition` exactly, so the
 * request body is the condition record and nothing is translated in between.
 */
export interface Conditions {
  /** K+, mM */
  k: number;
  /** Na+, mM */
  na: number;
  /** Li+ / NH4+, mM. Li+ is the standard G4-destabilising control. */
  li_nh4: number;
  /** Mg2+, mM */
  mg: number;
  ph: number;
  /** degrees C */
  temperature: number;
  /** % w/v PEG or equivalent */
  crowder_pct: number;
  /** strand concentration, uM */
  strand_conc: number;
}

/** What kind of statement a number is. Mirrors QuadCond's `target_semantics`. */
export type ClaimBasis =
  /** a physical measurement of the structure being predicted */
  | "measured"
  /** antibody occupancy at a genomic locus -- never P(folds) */
  | "genomic-proxy"
  /** no trained head covers this fold; condition-calibrated energy function */
  | "heuristic";

/**
 * Whether a value may be shown at all, and why not.
 *
 * `inDomain: false` does not mean "slightly uncertain". It means the query is
 * outside what the head was trained on, and the number beside it is an
 * extrapolation the interface hides until the user asks for it.
 */
export interface Applicability {
  inDomain: boolean;
  warnings: string[];
}

export interface Nucleotide {
  index: number;
  base: string;
  x: number;
  y: number;
  z: number;
  paired: boolean;
  role: "tetrad" | "loop" | "stem" | "ss" | "junction" | "pair";
  strand: number;
}

export interface EnsembleMember {
  id: string;
  kind: StructureKind;
  topology: string;
  probability: number;
  deltaG: number;
  predictedTm?: number;
  notation: string;
  details: string;
  nucleotides: Nucleotide[];
  strands: number[][];
  tetradBonds?: [number, number][];

  /**
   * Where this member's numbers came from. A fold whose dG is the local energy
   * function is `heuristic` and must be rendered differently from one whose Tm
   * came from a head trained on melting curves -- they are not the same kind of
   * statement, and the interface is the last place that distinction survives.
   */
  basis: ClaimBasis;
  /** The QuadCond head that supplied deltaG / predictedTm, when one did. */
  source?: string;
  applicability?: Applicability;
}

/**
 * Antibody occupancy at genomic loci: BG4 / iMab CUT&Tag.
 *
 * Deliberately NOT an ensemble member. These scores are an experimental
 * observation of a different quantity -- whether an antibody bound a site in
 * HEK293T cells -- and mixing them into a Boltzmann ensemble would state that
 * they are free energies of folding. They get their own panel.
 */
export interface GenomicEvidence {
  head: string;
  label: string;
  score: number;
  claim: string;
  applicability?: Applicability;
}

export interface ProteinBinder {
  id: string;
  name: string;
  gene: string;
  uniprot: string;
  kinds: StructureKind[];
  polymers: Polymer[];
  score: number;
  rationale: string;
  sequence: string;
  organism: string;
}

export interface AptamerTarget {
  id: string;
  name: string;
  description: string;
  preferredStructure: StructureKind;
  preferredPolymer: Polymer;
  targetKd: string;
}

export interface PredictionResult {
  polymer: Polymer;
  sequence: string;
  conditions: Conditions;
  ensemble: EnsembleMember[];
  proteins: ProteinBinder[];
  /** Rendered in its own panel, never weighted into the ensemble. */
  genomicEvidence?: GenomicEvidence[];
  /** True when the numbers came from the QuadCond service rather than locally. */
  fromService?: boolean;
  /**
   * Set when the tool declines to answer at all. An empty ensemble with a stated
   * reason, rather than a populated one with a warning: the viewer renders this
   * as a refusal, and there is no number on screen to misread.
   */
  refusal?: string;
}

/**
 * Nominal intracellular chemistry. Every field is a *choice*, not a measurement,
 * which is why QuadCond flags all eight as imputed when a caller supplies none.
 */
export const DEFAULT_CONDITIONS: Conditions = {
  k: 140,
  na: 10,
  li_nh4: 0,
  mg: 1,
  ph: 7.4,
  temperature: 37,
  crowder_pct: 0,
  strand_conc: 5,
};

/** Total monovalent cation, for the energy heuristics that only see ionic strength. */
export function monovalent(c: Conditions): number {
  return c.k + c.na + c.li_nh4;
}

export const KIND_LABEL: Record<StructureKind, string> = {
  duplex: "Duplex",
  hairpin: "Hairpin",
  "g-quadruplex": "G-quadruplex",
  "i-motif": "i-motif",
  "ac-motif": "AC-motif",
  cruciform: "Cruciform",
  triplex: "Triplex",
  "ss-coil": "Unfolded coil",
};
