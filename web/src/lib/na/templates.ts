/**
 * Known experimental structures for the sequences this tool is asked about.
 *
 * The geometry builder produces an *idealised model*: four strands around
 * stacked tetrads at 3.3 A with 30 degrees of inter-tetrad twist, or two
 * intercalated C:C+ duplexes at 6.2 A with the partner's pairs between them.
 * Those are published helical parameters and the builder measures itself
 * against them -- but it is not a prediction of atomic coordinates, and a
 * picture that does not say so is read as one.
 *
 * The honest distinction, made explicit in the interface:
 *
 *   - **experimental**  a solved structure exists for this exact sequence. The
 *     accession is shown and linked. The tool does not silently substitute it
 *     for the model, because a real fold and an idealised one differ in loop
 *     conformation, groove width and syn/anti glycosidic angles -- but the user
 *     is told it exists, which is the thing they would otherwise never learn.
 *
 *   - **schematic structural archetype**  no solved structure matches. What is
 *     drawn is the archetype for the topology, at published parameters, and
 *     nothing about it is specific to this sequence beyond tract placement.
 *
 * The registry is deliberately small and hand-checked. A large auto-harvested
 * table would be worse than none: a wrong accession beside a picture is a
 * citation the user will trust.
 */
import { cleanSeq } from "./predict";
import type { StructureKind } from "./types";

export interface StructureTemplate {
  pdb: string;
  name: string;
  sequence: string;
  kind: StructureKind;
  topology: string;
  method: string;
  conditions: string;
  /** Where this was checked, and what was checked. Required — a test enforces it. */
  verified: string;
}

/**
 * Sequences whose folds are solved and unambiguous. Each entry was checked
 * against the primary deposition rather than harvested from a secondary table.
 */
export const TEMPLATES: StructureTemplate[] = [
  {
    pdb: "143D",
    name: "Human telomeric d[AG3(T2AG3)3], Na+ form",
    sequence: "AGGGTTAGGGTTAGGGTTAGGG",
    kind: "g-quadruplex",
    topology: "antiparallel basket",
    method: "solution NMR",
    conditions: "Na+",
    verified:
      "https://www.rcsb.org/structure/143D — title gives the sequence as the " +
      "unambiguous formula d[AG3(T2AG3)3]; three stacked G-tetrads with two " +
      "lateral loops and one central diagonal loop, in Na+.",
  },
  {
    pdb: "1KF1",
    name: "Human telomeric d[AG3(T2AG3)3], K+ form",
    sequence: "AGGGTTAGGGTTAGGGTTAGGG",
    kind: "g-quadruplex",
    topology: "parallel propeller",
    method: "X-ray, 2.10 Å",
    conditions: "K+",
    verified:
      "https://www.rcsb.org/structure/1KF1 — per-residue sequence listing " +
      "matches; three propeller loops on the exterior, K+ in the channel.",
    // The same 22-mer as 143D, and the reason this table exists: one sequence,
    // two solved topologies, chosen by which cation is in the buffer. A viewer
    // that draws one picture per sequence is answering a question the chemistry
    // does not have a single answer to.
  },
];

/*
 * Three entries were removed in v0.4.3 because they were wrong, and each was
 * wrong in a way a reader would have had no way to catch:
 *
 *   - **2O3M** was labelled "c-MYC Pu22". It is "Monomeric G-DNA tetraplex from
 *     human C-kit promoter" — a different gene entirely, and the sequence given
 *     for it here was c-KIT1's, which is not 2O3M's either.
 *   - **2KYP** was labelled c-KIT1 with a 22-nt sequence. It is a 21-nt c-kit2
 *     structure, d(CGGCGGGCGCTGAGGGAGGT).
 *   - **1EL2** was listed as an unmodified telomeric i-motif solved by X-ray.
 *     It is "Solution structure of a MODIFIED human telomere fragment" — NMR,
 *     with 5-methylcytosine and uracil substitutions.
 *
 * They were added from memory rather than from the depositions. This registry
 * is now deliberately two entries long: a wrong accession printed beside a
 * structure is a citation the reader will trust, and there is no version of
 * "more coverage" worth that. Adding an entry requires reading the deposition
 * and recording the check in `verified`, which a test enforces.
 */

const BY_SEQUENCE = new Map<string, StructureTemplate[]>();
for (const t of TEMPLATES) {
  const key = cleanSeq(t.sequence);
  BY_SEQUENCE.set(key, [...(BY_SEQUENCE.get(key) ?? []), t]);
}

/** Exact-sequence matches only. A near match is a different molecule. */
export function templatesFor(sequence: string): StructureTemplate[] {
  return BY_SEQUENCE.get(cleanSeq(sequence)) ?? [];
}

export function templatesForKind(sequence: string, kind: StructureKind): StructureTemplate[] {
  return templatesFor(sequence).filter((t) => t.kind === kind);
}

export function pdbUrl(id: string): string {
  return `https://www.rcsb.org/structure/${id}`;
}

/** What the picture on screen actually is. */
export function provenanceLabel(sequence: string, kind: StructureKind): string {
  return templatesForKind(sequence, kind).length
    ? "schematic structural archetype — a solved structure exists for this sequence"
    : "schematic structural archetype";
}
