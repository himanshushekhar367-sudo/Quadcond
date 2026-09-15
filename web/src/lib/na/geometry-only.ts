/**
 * The local path: geometry, and nothing that looks like a prediction.
 *
 * `predict.ts` still exists and still contains a complete energy function --
 * hairpin, duplex, G4, i-motif, coil -- and that is precisely why nothing in
 * the running application may call it. Three things were true of the numbers it
 * produced, each fatal on its own:
 *
 *   1. **It could not tell potassium from sodium.** Every term reads total
 *      monovalent concentration, so 150 mM K+ and 150 mM Na+ returned a
 *      byte-identical ensemble -- while the measured difference on this very
 *      sequence is about 13 degC and the trained head reproduces it out-of-fold.
 *      The interface had six ion sliders in front of a function that saw one.
 *
 *   2. **The dG calibration was never applied here.** The span for Tel22 runs
 *      -27.9 to +18.1 kcal/mol: 46 kcal/mol, or 74 RT, which saturates the
 *      Boltzmann exponent clamp and pins the result at 0.5/0.5/0.000. The
 *      release notes describing a drop from 78 RT to 4.1 RT were describing the
 *      *service* path only.
 *
 *   3. **Its reference state was not zero.** The coil came out at +3.0 kcal/mol
 *      for a 22-mer, so every other free energy in the list floated on an offset
 *      of unstated size.
 *
 * Meanwhile the status bar said "geometry only -- nothing is substituted" and
 * the panel underneath rendered those probabilities as bars. This module is what
 * makes that sentence true: offline, the app builds coordinates and draws them,
 * and there is no number on screen for the absent service to be substituted for.
 */
import { geometryFor, type G4Topology } from "./geometry";
import { cleanSeq } from "./predict";
import { foldEnabled, RNA_SUPPORT } from "./capabilities";
import type { ClaimBasis, Conditions, Nucleotide, Polymer, StructureKind } from "./types";

/**
 * Which archetypes a strand can carry, by motif content alone.
 *
 * Motif detection is a regex over the sequence, not a thermodynamic claim, so it
 * is honest work for the client to do. Ranking them is not, and this returns
 * them in a fixed order rather than a scored one.
 */
function drawableKinds(seq: string): StructureKind[] {
  const kinds: StructureKind[] = [];
  if (/(?:G{3,}\w{1,12}){3,}G{3,}/i.test(seq)) kinds.push("g-quadruplex");
  if (/(?:C{3,}\w{1,12}){3,}C{3,}/i.test(seq)) kinds.push("i-motif");
  return kinds.filter(foldEnabled);
}

function revcomp(seq: string): string {
  const m: Record<string, string> = { A: "T", T: "A", G: "C", C: "G", N: "N" };
  return seq.split("").reverse().map((b) => m[b] ?? "N").join("");
}

/**
 * A drawable structure with no numbers attached.
 *
 * `deltaG` and `probability` are absent from the type by construction rather
 * than nullable: an optional number invites `?? 0` downstream, and a zero
 * renders as a bar.
 *
 * `strand` is not decoration. A G-quadruplex and an i-motif at one duplex
 * position sit on opposite strands, so a member is only identified by kind
 * *and* strand — and the previous version built members from the typed strand
 * alone, meaning a reverse-strand service card had no matching geometry and the
 * viewer silently fell back to `members[0]`. The card said "i-motif, strand −"
 * while the picture showed a parallel G4.
 */
export interface GeometryOnlyMember {
  id: string;
  kind: StructureKind;
  topology: string;
  strand: "+" | "-";
  strandSequence: string;
  notation: string;
  description: string;
  nucleotides: Nucleotide[];
  strands: number[][];
  tetradBonds?: [number, number][];
  basis: ClaimBasis;
}

/**
 * The geometry result.
 *
 * There is deliberately no `ensemble` field. It survived as `ensemble: []` for
 * a release after the local ensemble was removed, which kept every stale
 * `result.ensemble` read compiling — `ViewerPane` blanked the 3D pane and
 * `Controls` lost structure provenance and dot-bracket notation, both without a
 * type error. An empty array is a valid array; an absent field is not.
 */
export interface GeometryResult {
  polymer: Polymer;
  sequence: string;
  complement: string;
  conditions: Conditions;
  members: GeometryOnlyMember[];
  refusal?: string;
}

/** DNA unless the sequence contains U, or the caller insists. */
function polymerOf(seq: string, hint: Polymer): Polymer {
  if (/U/i.test(seq)) return "RNA";
  return hint;
}

/**
 * Dot-bracket for the drawn archetype. Pure string work over the motif regex —
 * it describes the picture and asserts nothing about stability.
 */
function notationFor(kind: StructureKind, seq: string): string {
  if (kind === "g-quadruplex")
    return seq.replace(/G+/g, (m) => "(".repeat(m.length)).replace(/[^G()]/g, ".");
  if (kind === "i-motif")
    return seq.replace(/C+/g, (m) => "[".repeat(m.length)).replace(/[^C[\]]/g, ".");
  if (kind === "duplex") return "(".repeat(seq.length);
  return ".".repeat(seq.length);
}

const DESCRIPTION =
  "Schematic archetype at published helical parameters. Which structure a " +
  "sequence actually adopts, and how stable it is, come from QuadCond — they " +
  "are not computed here.";

export function geometryOnly(
  raw: string,
  polymerHint: Polymer,
  conditions: Conditions,
): GeometryResult {
  const sequence = cleanSeq(raw);
  const polymer = polymerOf(sequence, polymerHint);
  const complement = revcomp(sequence);
  const base = { polymer, sequence, complement, conditions, members: [] as GeometryOnlyMember[] };

  if (polymer === "RNA" && !RNA_SUPPORT.enabled) {
    return { ...base, refusal: RNA_SUPPORT.reason };
  }
  if (sequence.length < 6) {
    return { ...base, refusal: "Sequence too short to build a structure." };
  }

  const members: GeometryOnlyMember[] = [];

  // Both strands, because the service answers about both. A G-rich input has
  // its i-motif on the complement, and a picture must exist for that card.
  for (const [strand, strandSequence] of [
    ["+", sequence],
    ["-", complement],
  ] as const) {
    for (const kind of drawableKinds(strandSequence)) {
      // All three topologies, including hybrid (3+1). Emitting only two meant a
      // hybrid card — which is what the heads predict for human telomeric DNA
      // in K+, the single most common query this tool gets — fell through to
      // the parallel archetype.
      const topologies: (G4Topology | undefined)[] =
        kind === "g-quadruplex" ? ["parallel", "antiparallel", "hybrid"] : [undefined];
      for (const topo of topologies) {
        const geo = geometryFor(kind, strandSequence, topo);
        members.push({
          id: `${kind}${topo ? `-${topo}` : ""}@${strand}`,
          kind,
          topology: topo ?? kind,
          strand,
          strandSequence,
          notation: notationFor(kind, strandSequence),
          description: DESCRIPTION,
          nucleotides: geo.nucleotides,
          strands: geo.strands,
          tetradBonds: geo.tetradBonds,
          basis: "heuristic",
        });
      }
    }
  }

  // The competing states, drawn once for the duplex position as a whole.
  for (const kind of (["duplex", "ss-coil"] as StructureKind[]).filter(foldEnabled)) {
    const geo = geometryFor(kind, sequence, undefined);
    members.push({
      id: `${kind}@+`,
      kind,
      topology: kind,
      strand: "+",
      strandSequence: sequence,
      notation: notationFor(kind, sequence),
      description: DESCRIPTION,
      nucleotides: geo.nucleotides,
      strands: geo.strands,
      tetradBonds: geo.tetradBonds,
      basis: "heuristic",
    });
  }

  return { ...base, members };
}

/**
 * Find the drawable member for a service card.
 *
 * Matched on kind **and strand**, then topology. Returns undefined rather than
 * a fallback: showing the wrong molecule is worse than showing none, and a
 * silent `members[0]` is how the card and the picture came to disagree.
 */
export function memberForCard(
  members: GeometryOnlyMember[],
  card: { kind: string; strand?: "+" | "-"; topology?: string },
): GeometryOnlyMember | undefined {
  const sameKind = members.filter((m) => m.kind === card.kind);
  const pool = card.strand ? sameKind.filter((m) => m.strand === card.strand) : sameKind;
  if (!pool.length) return undefined;

  // Whole-word match, not substring: "antiparallel basket".includes("parallel")
  // is true, so a substring test picked the parallel archetype for an
  // antiparallel card — and picked it first, because parallel is emitted first.
  const words = new Set((card.topology ?? "").toLowerCase().split(/[^a-z]+/).filter(Boolean));
  const named = pool.find((m) => words.has(m.topology.toLowerCase()));
  if (named) return named;

  // Nothing matched. If this kind has more than one archetype, there is a
  // choice to be made and no basis for making it — so make none. The `??
  // pool[0]` that used to sit here is how a hybrid card came to be drawn as a
  // parallel quadruplex: parallel is emitted first, so "no match" silently
  // meant "parallel". A named topology is a specific claim about a specific
  // fold, and substituting a different one answers a question nobody asked.
  if (pool.length > 1) return undefined;

  // Exactly one archetype for this kind (an i-motif card, say) — no ambiguity.
  return pool[0];
}
