import type { Polymer } from "./types";

export interface ExampleSeq {
  /**
   * What this tool can actually say about the preset.
   *
   * A preset list is a set of promises. This one used to advertise "FMRP /
   * DHX36 binders" for a sequence the service rejects outright, and Mg2+
   * occupancy for an AC-motif no head models. A single `unsupported` boolean
   * was not enough either: the AC-motif oligo carries four C-tracts, so its
   * i-motif card is real evidence and only the AC interpretation is missing —
   * badging the whole thing "not covered" is its own kind of wrong.
   *
   *   full     — the structures drawn are the ones the heads predict
   *   partial  — some structure here is covered, the named one is not
   *   none     — nothing relevant is modelled
   *   refused  — the service will not answer at all
   */
  coverage?: "full" | "partial" | "none" | "refused";
  id: string;
  name: string;
  polymer: Polymer;
  sequence: string;
  note: string;
}

export const EXAMPLES: ExampleSeq[] = [
  {
    id: "telo-dna",
    name: "Human telomere",
    polymer: "DNA",
    sequence: "GGGTTAGGGTTAGGGTTAGGG",
    note: "Parallel / hybrid G4 in K+; C-rich complement forms i-motif at low pH.",
  },
  {
    id: "cmyc",
    name: "c-MYC Pu27",
    polymer: "DNA",
    sequence: "TGGGGAGGGTGGGGAGGGTGGGGAAGG",
    note: "Oncogene promoter G4. High K+ locks a parallel fold.",
  },
  {
    id: "telo-c",
    name: "Telomeric C-strand",
    polymer: "DNA",
    sequence: "CCCTAACCCTAACCCTAACCCT",
    note: "Classic intramolecular i-motif. Occupancy collapses above pH ~6.5.",
  },
  {
    id: "terra",
    name: "TERRA RNA",
    polymer: "RNA",
    sequence: "GGGUUAGGGUUAGGGUUAGGG",
    note: "Telomeric RNA. The atlas is entirely DNA, so the service refuses it. Kept so the refusal is visible rather than hypothetical.",
    coverage: "refused",
  },
  {
    id: "hairpin",
    name: "GC hairpin",
    polymer: "DNA",
    sequence: "GCGCGCGCAAAAGCGCGCGC",
    note: "Stem-loop. No head predicts hairpin stability and there is no hairpin archetype, so only duplex and coil are drawn.",
    coverage: "none",
  },
  {
    id: "acmotif",
    name: "AC-motif oligo",
    polymer: "DNA",
    sequence: "AACCCAACCCAACCCAACCCAA",
    note: "Four C-tracts, so the i-motif card is real. The A+/C intercalation this construct is named for is NOT modelled — no head covers it and no atlas row constrains it.",
    coverage: "partial",
  },
  {
    id: "cruciform",
    name: "Palindromic cruciform",
    polymer: "DNA",
    sequence: "CGCGAATTCGCGTTTTCGCGAATTCGCG",
    note: "Inverted repeats. Cruciform extrusion needs superhelical density, for which this tool has no term. Only duplex and coil are drawn.",
    coverage: "none",
  },
  {
    id: "mixed",
    name: "Competing G/C oligo",
    polymer: "DNA",
    sequence: "GGGTTAGGGCCCTAACCCGGGTTAGGG",
    note: "Carries both a G-tract and a C-tract system, so both strand cards appear. Side by side is not a competition — the duplex is not modelled.",
  },
  {
    id: "rna-hp",
    name: "RNA tetraloop",
    polymer: "RNA",
    sequence: "GGGCGCAAGCGCC",
    note: "RNA hairpin. No rG4 or RNA measurement stands behind any head.",
    coverage: "refused",
  },
];
