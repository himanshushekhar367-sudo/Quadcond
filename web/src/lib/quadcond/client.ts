/**
 * Client for the QuadCond prediction service.
 *
 * The viewer used to compute its own numbers from a hand-written energy
 * function. Those numbers had no referent: the free energies were roughly ten
 * times too steep relative to RT, so the Boltzmann step had no dynamic range,
 * and benchmarked against 2,274 measured melting temperatures the heuristic
 * scored R² −0.465 — worse than predicting the mean.
 *
 * So the app no longer computes what a trained model can. It asks — and when the
 * service cannot be reached it shows no numbers at all rather than falling back
 * to the local function, which could not tell potassium from sodium and returned
 * a byte-identical answer for both.
 *
 * The service has no accounts, no sessions and no cookies, which is a
 * requirement rather than a convenience: a web server that makes users register
 * is not eligible for the NAR Web Server Issue, and a prediction endpoint has
 * nothing to authenticate anyway.
 */

/**
 * The provenance block every comparison workflow returns.
 *
 * One declared shape rather than `Record<string, unknown>` per endpoint. Until
 * v0.4.8 only `/batch` carried the model artifact hash and the atlas
 * fingerprint, so a mutation or condition export could name a model *version*
 * and not the file that produced it — which is the difference between a figure
 * that can be regenerated and one that cannot.
 */
export type RunRecord = {
  workflow: "mutation_scan" | "condition_scan" | "batch";
  quadcond_model_version: string | null;
  model_artifact_sha256: string | null;
  atlas_fingerprint: string | null;
  prediction_schema_version: string;
  heads_requested?: string[] | "all";
  strand_by_head?: Record<string, "+" | "-">;
  [key: string]: unknown;
};

export type ClaimBasis = {
  has_experimental_observation?: boolean;
  target_semantics: string;
  biophysically_grounded: boolean;
  calibration_scope?: string;
};

export type Applicability = {
  in_domain: boolean;
  warnings: string[];
  biophysically_grounded?: boolean;
};

export type EvidenceCard = {
  id: string;
  kind: string;
  topology: string;
  strand: "+" | "-";
  strand_sequence: string;
  source: string;
  predictedTm?: number;
  predictedTm_interval?: [number, number];
  predictedPhT?: number;
  predictedPhT_interval?: [number, number];
  foldedFraction?: number;
  /** Why the folded fraction is not an occupancy. Rendered, not dropped. */
  foldedFraction_note?: string;
  // No `probabilityOfFolding`. The field is gone from the service — the
  // discrimination scores it used to carry now sit in their own block with
  // their own provenance — and leaving the optional declaration behind is the
  // `ensemble: []` mistake again: a type that still admits the field is a type
  // a component can be written against, and the next thing to render it will
  // compile. A branch of this project did exactly that, and went further,
  // normalising card scores plus a hardcoded 0.40 hairpin and 0.95 duplex into
  // "population percentages" that summed to one. Deleting the field is what
  // makes that a compile error rather than a code review.
  claim_basis: ClaimBasis;
  applicability: Applicability;
  motifs?: number;
};

export type RefusedReadout = {
  refused: true;
  refusal_reason: string;
};

export type GenomicEvidenceItem = {
  head: string;
  label: string;
  strand: "+" | "-";
  strand_sequence: string;
  claim: string;
  applicability: Applicability;
} & (RefusedReadout | { refused?: false; score: number });

export type LocusOverlap = {
  applicability: Applicability;
  strand: "+" | "-";
  strand_sequence: string;
  note: string;
} & (RefusedReadout | {
  refused?: false;
  posterior: Record<string, number>;
  argmax: string;
  claim_basis: ClaimBasis;
});

export type EvidenceResponse = {
  quadcond_version: string;
  model_version: string;
  prediction_schema_version?: string;
  sequence: string;
  complement: string;
  condition: string;
  condition_detail: Record<string, number>;
  condition_imputed_fields: string[];
  /**
   * Absent since v0.4.3 — the Boltzmann step it scaled was withdrawn. Declared
   * optional rather than deleted so an older service is still readable, and
   * never required, so no component can be written against a field the current
   * service does not send.
   */
  RT_kcal?: number;
  evidence: EvidenceCard[];
  /**
   * Antibody occupancy scores. Rendered in their own panel and never weighted
   * into the ensemble: putting a CUT&Tag score in a partition function states
   * that it is a folding free energy, which is the exact misreading the whole
   * claim-basis apparatus exists to prevent.
   */
  genomic_evidence?: GenomicEvidenceItem[];
  genomic_evidence_note?: string;
  locus_peak_overlap_state: LocusOverlap | null;
  evidence_note: string;
  /**
   * Present when the service is serving a model trained on generated data.
   *
   * `/evidence` was the one numerical endpoint that did not spread the banner
   * into its response, so the same predictor produced a flagged `/predict`
   * payload and an unflagged `/evidence` payload carrying the same numbers in
   * cards — and `/evidence` is what the viewer renders by default. "Every
   * response is flagged" has to include this one.
   */
  synthetic_model?: SyntheticBanner;
  model_artifact_sha256?: string | null;
  atlas_fingerprint?: string | null;
  /** Also absent since v0.4.3: shipping it would let a caller rebuild the
   *  partition function that was withdrawn. */
  calibration?: Record<string, unknown>;
};

/**
 * What one head can and cannot respond to, as the backend reports it.
 *
 * The four states are not decoration. `varied` means training data spans the
 * axis and a prediction can depend on it. `fixed` means every training record
 * used one value, so the head has no response to apply — an absence, not an
 * extrapolation. `predicts` means the head outputs that axis rather than
 * consuming it (a Tm head does not take a temperature). `unknown` means no
 * range was recorded, treated as `fixed` because assuming a response nobody
 * measured is the failure this exists to prevent.
 *
 * The viewer used to show six ion sliders in front of heads that consume one or
 * none of them. A slider that moves while a number sits still teaches the reader
 * that the head is *insensitive* to that variable, which is a scientific claim
 * nobody made.
 */
export type AxisState = {
  state: "varied" | "fixed" | "predicts" | "unknown";
  label: string;
  note: string;
  value?: number;
  min?: number;
  max?: number;
  n_unique?: number;
};

export type HeadCapabilities = {
  workspace: "sequence_evidence" | "condition_response" | "genomic_context" | "diagnostics";
  workspace_label: string;
  axes: Record<string, AxisState>;
  responds_to: string[];
  inert_axes: string[];
  predicted_axes: string[];
  use_conditions: boolean;
};

export type HeadInfo = ClaimBasis &
  HeadCapabilities & {
    task: string;
    kind: string;
    target: string;
    units: string;
    claim: string;
    n_training_rows?: number;
    metrics?: Record<string, unknown>;
    applicability?: Record<string, unknown>;
    target_semantics_label?: string;
  };

export type Workspace = {
  id: HeadCapabilities["workspace"];
  label: string;
  note: string;
  heads: string[];
};

export type ServiceInfo = {
  quadcond_version: string;
  model_version: string;
  headline: string;
  heads: Record<string, HeadInfo>;
  /** The four groups, served populated. An empty one is still sent. */
  workspaces: Workspace[];
  condition_axes: { id: string; label: string }[];
  presets?: Record<string, Record<string, number>>;
  atlas_composition?: Record<string, unknown> | null;
  prediction_schema_version?: string;
  synthetic_model?: SyntheticBanner;
};

/**
 * Present only when the service is serving a model trained on generated data.
 *
 * It rides on every response rather than being printed at startup, because a
 * JSON file that ends up feeding a figure has to carry its own warning.
 */
export type SyntheticBanner = {
  provenance: "synthetic";
  note: string;
  warning: string;
};

export type ServiceHealth = {
  status: "ok" | "unavailable";
  quadcond_version?: string;
  assets?: { name: string; state: string; path: string | null }[];
  authentication?: string;
  error?: string;
};

export type QuadcondConditions = {
  k?: number;
  na?: number;
  li_nh4?: number;
  mg?: number;
  ph?: number;
  temperature?: number;
  crowder_pct?: number;
  strand_conc?: number;
  preset?: string;
};

const DEFAULT_BASE =
  (typeof import.meta !== "undefined" &&
    (import.meta as { env?: Record<string, string> }).env?.VITE_QUADCOND_URL) ||
  "http://127.0.0.1:8765";

export class QuadcondUnavailable extends Error {
  constructor(message: string) {
    super(message);
    this.name = "QuadcondUnavailable";
  }
}

async function post<T>(base: string, path: string, body: unknown,
                       signal?: AbortSignal): Promise<T> {
  let res: Response;
  try {
    res = await fetch(`${base}${path}`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      // no credentials: there is no session, and sending them would be the
      // first step towards a server that needs one
      credentials: "omit",
      body: JSON.stringify(body),
      signal,
    });
  } catch (cause) {
    throw new QuadcondUnavailable(
      `cannot reach the QuadCond service at ${base}. Start it with ` +
        `'python -m quadcond.service', or set VITE_QUADCOND_URL.`,
    );
  }
  const payload = (await res.json().catch(() => ({}))) as Record<string, unknown>;
  if (!res.ok) {
    throw new QuadcondUnavailable(
      typeof payload.error === "string" ? payload.error : `service returned ${res.status}`,
    );
  }
  return payload as T;
}

export async function health(base = DEFAULT_BASE): Promise<ServiceHealth> {
  try {
    const res = await fetch(`${base}/health`, { credentials: "omit" });
    return (await res.json()) as ServiceHealth;
  } catch {
    return { status: "unavailable", error: `no service at ${base}` };
  }
}

export async function info(base = DEFAULT_BASE) {
  const res = await fetch(`${base}/info`, { credentials: "omit" });
  if (!res.ok) throw new QuadcondUnavailable(`service returned ${res.status}`);
  return res.json();
}

export async function evidence(
  sequence: string,
  conditions: QuadcondConditions = {},
  base = DEFAULT_BASE,
  signal?: AbortSignal,
): Promise<EvidenceResponse> {
  return post<EvidenceResponse>(base, "/evidence", { sequence, ...conditions }, signal);
}

export async function predict(
  sequences: string[],
  conditions: QuadcondConditions = {},
  base = DEFAULT_BASE,
  signal?: AbortSignal,
) {
  return post(base, "/predict", { sequences, ...conditions, neighbours: 3 }, signal);
}

// --------------------------------------------------------------- comparisons

/** One cell of a mutation scan: what one head says about one substitution. */
export type MutationCell = {
  /** Which strand this head was asked about. A G4 and an i-motif at one duplex
   *  position sit on opposite strands, so without this a reader comparing two
   *  rows is comparing two sequences. */
  strand: "+" | "-";
  motif: {
    state:
      | "motif_retained"
      | "motif_lost"
      | "motif_gained"
      | "motif_count_changed"
      | "no_motif";
    before: number;
    after: number;
    note?: string;
  };
  wild_type: { state: string; reason?: string; warnings?: string[] };
  mutant: { state: string; reason?: string; warnings?: string[] };
  /** Present only when both sides carry a value and the head is applicable to
   *  both. There is deliberately no key to read a zero out of. */
  wild_type_value?: number;
  mutant_value?: number;
  delta?: number;
  /** Spread of the paired per-estimator differences — not the marginal
   *  intervals propagated, and not a validated interval on the difference. */
  delta_paired_sd?: number;
  n_estimators?: number;
  // No `delta_note` here. It is one sentence for the whole scan and it lives at
  // the top level of the response — 330 copies of it in a 22-mer's payload were
  // most of the bytes and no more informative for being repeated. Declaring it
  // per cell meant the panel read `rows[0].heads[h].delta_note`, found
  // undefined, and silently fell back to a shorter warning that omitted the
  // actual limitation. A type that still admits a field the service does not
  // send is a type a component can be written against.
};

export type MutationRow = {
  position: number;
  position_1based: number;
  wild_type_base: string;
  mutant_base: string;
  label: string;
  sequence: string;
  complement: string;
  heads: Record<string, MutationCell>;
};

export type MutationScan = {
  quadcond_version: string;
  model_version: string;
  scan: "mutation";
  /** Frozen at the top of the scan and returned with it, so a stored table can
   *  never be paired with a different baseline without the mismatch showing. */
  wild_type: {
    sequence: string;
    complement: string;
    length: number;
    condition: string;
    condition_detail: Record<string, number>;
    condition_imputed_fields: string[];
    predictions: Record<string, Record<string, unknown> | null>;
    strand_by_head: Record<string, "+" | "-">;
    /** The two strand strings once, rather than on every cell. */
    strand_sequence: { "+": string; "-": string };
  };
  heads: Record<string, HeadInfo & { units: string }>;
  substitutions: MutationRow[];
  /** Said once, for the whole scan: what the paired spread is and is not. */
  delta_note: string;
  ranking_note: string;
  run_record: RunRecord;
  atlas_fingerprint?: string | null;
  schema_version?: string;
  synthetic_model?: SyntheticBanner;
};

export type ConditionPoint = {
  value_of_axis: number;
  /** Which strand this point was predicted on. Two series in one panel can
   *  otherwise be about two different molecules. */
  strand: "+" | "-";
  applicability: { state: string; reason?: string; warnings?: string[] };
  value?: number;
  interval?: [number, number];
  probability?: number;
  folded_fraction?: number;
  transition_model?: string;
};

/**
 * Did the prediction actually move?
 *
 * `responds_to` says an axis varied in training, which licenses a dependence
 * but does not demonstrate one — a composition-matched shuffle inherits its
 * parent row's buffer, so a classification head shows wide training variation
 * in an axis its label never depended on. This is the measurement, taken
 * against the head's own error scale, because a chart auto-scales and a 0.2 °C
 * wobble fills the panel exactly like a 20 °C swing.
 */
export type ObservedResponse = {
  /**
   * `below_error_scale` replaces the old `flat`. "Flat" reads as a measured
   * insensitivity to the axis, which is a claim nobody made; what was actually
   * computed is a response range set beside a marginal residual half-width,
   * which is descriptive and is not a test that a difference is unresolvable.
   */
  verdict:
    | "responds"
    | "below_error_scale"
    | "unknown_scale"
    | "refused_throughout"
    | "insufficient_points"
    | "insufficient_supported_points";
  quantity?: "value" | "probability";
  min?: number;
  max?: number;
  spread?: number;
  error_scale?: number | null;
  /** Always `"in_domain_points_only"`. Out-of-domain points are returned with
   *  their state but never enter the summary — an unsupported tail of a sweep
   *  deciding the headline verdict is the opposite of what a domain check is
   *  for. */
  summary_basis?: "in_domain_points_only";
  n_points?: number;
  n_supported?: number;
  n_out_of_domain?: number;
  n_refused?: number;
  note?: string;
};

export type ConditionSeries = {
  head: string;
  kind: string;
  task: string;
  target: string;
  units: string;
  claim: string;
  target_semantics: string;
  biophysically_grounded: boolean;
  n_training_rows?: number;
  training_span?: { min?: number; max?: number; n_unique?: number };
  /** The strand this head was asked about, and that strand's sequence. */
  strand: "+" | "-";
  strand_sequence: string;
  observed_response: ObservedResponse;
  points: ConditionPoint[];
};

export type ConditionScan = {
  quadcond_version: string;
  model_version: string;
  scan: "condition";
  sequence: string;
  complement: string;
  /** The two strand strings once, keyed by strand. */
  strand_sequence: { "+": string; "-": string };
  axis: string;
  axis_label: string;
  /** The axis values, and only those. This key used to arrive holding one
   *  head's list of point objects, because a loop inside the service reused the
   *  name that held it — so a client declaring `number[]` received `object[]`. */
  values: number[];
  base_condition: Record<string, number>;
  condition_imputed_fields?: string[];
  run_record: RunRecord;
  summary_note?: string;
  neighbour_tiers?: string[];
  series: Record<string, ConditionSeries>;
  /** Named with a reason rather than omitted: a head missing from a chart reads
   *  as "not applicable", which is one of two very different things. */
  inert_heads: { head: string; axis_state: string; reason: string; claim: string }[];
  measured_neighbours: Record<string, unknown>[];
  neighbour_note: string;
  synthetic_model?: SyntheticBanner;
};

export type BatchResult = {
  quadcond_version: string;
  scan: "batch";
  /** Counts are by input index, never by identifier: two records sharing a
   *  name are two sequences, and reporting them as one reads as a silent
   *  exclusion that never happened. */
  run_record: RunRecord & {
    input_mode?: "fasta" | "lines" | "records";
    n_sequences_in: number;
    n_sequences_scored: number;
    n_sequences_excluded: number;
    n_rows: number;
    duplicate_identifiers?: string[];
    reconciled?: boolean;
  };
  results: (Record<string, unknown> & { id: string; record_index?: number })[];
  /** Never silently dropped. A batch of 200 that returns 197 rows with no
   *  explanation is a batch whose three missing sequences get noticed after the
   *  figure is drawn, if at all. */
  excluded: { id: string; record_index?: number; reason: string; input?: string }[];
  synthetic_model?: SyntheticBanner;
};

/**
 * A substitution seen from two directions at once.
 *
 * The structural axis is one QuadCond head's delta for this edit; the
 * regulatory axis is whatever AlphaGenome Atlas says about the same edit at the
 * same coordinate. There is deliberately no field combining them: the two are
 * model output on different scales and no weighting between them has been
 * calibrated, so the response carries both, the quadrant they fall in against
 * stated thresholds, and a rank that is high only when both are high.
 */
export type VariantRow = {
  variant: { chromosome: string; position: number; reference: string; alternate: string };
  label: string;
  mutation: string;
  position_in_window: number;
  structural: {
    head: string | null;
    strand: "+" | "-" | null;
    delta?: number | null;
    delta_paired_sd?: number | null;
    motif_state: string | null;
    magnitude: number | null;
    rank?: number | null;
    applicability?: { state: string; warnings?: string[] } | null;
  };
  /** Null means Atlas has no record for this substitution — which is not a low
   *  score, and is why such a row has no quadrant. */
  regulatory: {
    scorer: string;
    score: number;
    rank?: number | null;
    all_scores: Record<string, number>;
    source: string;
  } | null;
  quadrant: "both" | "structural_only" | "regulatory_only" | "neither" | "unclassified";
  quadrant_reason?: string;
  combined_rank: number | null;
  heads: Record<string, MutationCell>;
};

/** A structural finding with no number by construction. Kept off the ranked
 *  table and out of reach of anything that would give it an invented delta. */
export type StructuralCandidate = {
  label: string;
  mutation: string;
  strand: "+" | "-" | null;
  motif_state: "motif_lost" | "motif_gained" | "motif_count_changed";
  regulatory_score: number | null;
  regulatory_rank: number | null;
  why_no_delta: string;
};

/**
 * Which model file answered, by its contents.
 *
 * `matches_manifest` is three-valued deliberately. `null` means the comparison
 * could not be made, and an unknown must not render as a pass.
 */
export type ModelIdentity = {
  model_version_reported_by_artifact: string | null;
  model_file: string | null;
  model_file_sha256: string | null;
  dataset_fingerprint: string | null;
  manifest_release: string | null;
  manifest_expected_model_sha256: string | null;
  matches_manifest: boolean | null;
  note: string;
};

/** Per-head, per-axis: did this head learn a response to this condition at all? */
export type ConditionResponsiveness = Record<
  string,
  Record<string, "varied" | "fixed" | "predicts" | "unknown">
>;

export type VariantScan = {
  quadcond_version: string;
  model_version: string;
  scan: "variant";
  run_record: RunRecord & {
    locus: { chromosome: string; start: number; end: number };
    regulatory_source: { source: string; note?: string; path?: string; scorers?: string[] };
    regulatory_error?: string | null;
    structural_axis_head?: string | null;
    regulatory_axis_scorer?: string | null;
    structural_threshold?: number | null;
    regulatory_threshold?: number | null;
    threshold_basis?: string;
    provenance?: { model: ModelIdentity; environment?: Record<string, unknown> };
    condition_responsiveness?: ConditionResponsiveness;
    condition_note?: string;
  };
  window: {
    chromosome: string; start: number; end: number; sequence: string;
    g4_elements: Record<string, unknown>[];
    im_elements: Record<string, unknown>[];
  };
  wild_type: MutationScan["wild_type"];
  heads: Record<string, HeadInfo & { units: string }>;
  variants: VariantRow[];
  quadrant_counts: Record<string, number>;
  structural_candidates: Record<string, StructuralCandidate[]>;
  structural_candidate_counts: Record<string, number>;
  candidate_note: string;
  quadrant_note: string;
  rank_note: string;
  delta_note: string;
  ranking_note: string;
  interpretation_note: string;
  synthetic_model?: SyntheticBanner;
};

export async function scanVariant(
  sequence: string,
  locus: { chromosome: string; start: number },
  conditions: QuadcondConditions = {},
  opts: {
    heads?: string[];
    structural_head?: string;
    regulatory_scorer?: string;
    structural_threshold?: number;
    regulatory_threshold?: number;
  } = {},
  base = DEFAULT_BASE,
  signal?: AbortSignal,
): Promise<VariantScan> {
  // The regulatory source is not a request parameter. The operator configures
  // one when the service starts, so a client cannot name a file path or post
  // an API key.
  return post<VariantScan>(base, "/scan/variant",
    { sequence, ...locus, ...conditions, ...opts }, signal);
}

export async function scanMutations(
  sequence: string,
  conditions: QuadcondConditions = {},
  opts: { heads?: string[]; positions?: number[] } = {},
  base = DEFAULT_BASE,
  signal?: AbortSignal,
): Promise<MutationScan> {
  return post<MutationScan>(base, "/scan/mutations",
    { sequence, ...conditions, ...opts }, signal);
}

export async function scanConditions(
  sequence: string,
  axis: string,
  span: { values?: number[]; min?: number; max?: number; points?: number },
  conditions: QuadcondConditions = {},
  base = DEFAULT_BASE,
  signal?: AbortSignal,
): Promise<ConditionScan> {
  return post<ConditionScan>(base, "/scan/conditions",
    { sequence, axis, ...span, ...conditions }, signal);
}

export async function batch(
  input: {
    fasta?: string;
    records?: { id: string; sequence: string }[];
    /** `auto` reads the text as FASTA when any line starts with ">", and as one
     *  sequence per line otherwise — which is what the input label promises. */
    parse_mode?: "auto" | "fasta" | "lines";
  },
  conditions: QuadcondConditions = {},
  opts: { heads?: string[]; conditions?: Record<string, number>[] } = {},
  base = DEFAULT_BASE,
  signal?: AbortSignal,
): Promise<BatchResult> {
  return post<BatchResult>(base, "/batch", { ...input, ...conditions, ...opts }, signal);
}

/**
 * Is this fold's number safe to present as biophysics?
 *
 * Used by the viewer to decide whether a value gets a plain label or a proxy
 * badge. The distinction is not cosmetic: a genomic-proxy score and a measured
 * melting temperature rendered identically next to a 3D model is how a tool
 * teaches somebody something untrue.
 */
export function isBiophysical(fold: EvidenceCard): boolean {
  return fold.claim_basis?.biophysically_grounded === true;
}

export function claimBadge(fold: EvidenceCard): { label: string; tone: "solid" | "warn" | "muted" } {
  const ts = fold.claim_basis?.target_semantics;
  // Not "measured". The number on the card is a model's prediction whose
  // training *labels* were measurements — which is a claim about the corpus,
  // not about this value. A badge reading "measured" beside a predicted melting
  // temperature is the single most misreadable word in this interface, and raw
  // measurements (the atlas neighbours) need a label of their own that this one
  // was occupying.
  if (ts === "biophysical") {
    return { label: "prediction · measured training data", tone: "solid" };
  }
  if (ts === "genomic_proxy") return { label: "genomic proxy", tone: "warn" };
  if (ts === "heuristic") return { label: "heuristic", tone: "warn" };
  if (ts === "reference") return { label: "reference", tone: "muted" };
  return { label: ts ?? "unknown", tone: "muted" };
}


/**
 * Folds the interface will show without being asked.
 *
 * A fold whose head reported `in_domain: false` is withheld by default. The
 * reference state is always kept: without it the remaining probabilities are
 * conditional on folding, which is a different question from the one the panel
 * appears to be answering.
 */
export function visibleFolds(
  folds: EvidenceCard[],
  showExtrapolation: boolean,
): EvidenceCard[] {
  if (showExtrapolation) return folds;
  return folds.filter((f) => f.applicability?.in_domain !== false || f.id === "unfolded");
}

export function extrapolatedFolds(folds: EvidenceCard[]): EvidenceCard[] {
  return folds.filter((f) => f.applicability?.in_domain === false && f.id !== "unfolded");
}

/** Every distinct applicability warning across an ensemble, deduplicated. */
export function allWarnings(folds: EvidenceCard[]): string[] {
  return [...new Set(folds.flatMap((f) => f.applicability?.warnings ?? []))];
}
