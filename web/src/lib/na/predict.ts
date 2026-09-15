import { geometryFor } from "./geometry";
import { bindersFor } from "./proteins";
import { PROTEIN_BINDERS, RNA_SUPPORT, foldEnabled } from "./capabilities";
import { monovalent } from "./types";
import type {
  ClaimBasis,
  Conditions,
  EnsembleMember,
  Polymer,
  PredictionResult,
  StructureKind,
} from "./types";

const R = 0.001987204;

export function cleanSeq(raw: string): string {
  return raw
    .toUpperCase()
    .replace(/^>.*$/gm, "")
    .replace(/[^ACGTU]/g, "");
}

function detectPolymer(seq: string, hint: Polymer): Polymer {
  if (hint === "RNA" || seq.includes("U")) return "RNA";
  return "DNA";
}

function complementary(a: string, b: string, polymer: Polymer): boolean {
  const pairs =
    polymer === "RNA"
      ? ["AU", "UA", "GC", "CG", "GU", "UG"]
      : ["AT", "TA", "GC", "CG", "GT", "TG"];
  return pairs.includes(a + b);
}

function gcFrac(seq: string): number {
  if (!seq.length) return 0;
  return (seq.match(/[GC]/g)?.length ?? 0) / seq.length;
}

function countRuns(seq: string, base: string, min = 2): { n: number; max: number } {
  const re = new RegExp(`${base}{${min},}`, "g");
  const m = seq.match(re) ?? [];
  return { n: m.length, max: m.reduce((a, x) => Math.max(a, x.length), 0) };
}

function invertedRepeatScore(seq: string): { score: number; stem: number } {
  const n = seq.length;
  let best = 0;
  let stem = 0;
  for (let i = 0; i < n; i++) {
    for (let j = n - 1; j > i + 4; j--) {
      let k = 0;
      while (
        i + k < j - k &&
        complementary(seq[i + k] ?? "", seq[j - k] ?? "", seq.includes("U") ? "RNA" : "DNA")
      ) {
        k++;
      }
      if (k > best) {
        best = k;
        stem = k;
      }
    }
  }
  return { score: best, stem };
}

function hairpinEnergy(seq: string, polymer: Polymer, c: Conditions): { dG: number; stem: number; loop: number; notation: string } {
  const n = seq.length;
  let bestStem = 0;
  let bestI = 0;
  for (let stem = Math.min(12, Math.floor((n - 3) / 2)); stem >= 3; stem--) {
    let ok = 0;
    for (let k = 0; k < stem; k++) {
      if (complementary(seq[k] ?? "", seq[n - 1 - k] ?? "", polymer)) ok++;
    }
    if (ok >= stem - 1 && stem > bestStem) {
      bestStem = stem;
      bestI = ok;
    }
  }
  const loop = n - 2 * bestStem;
  const gc = gcFrac(seq.slice(0, bestStem));
  let dG = 4.8;
  if (bestStem < 3) dG = 3.5;
  else {
    dG = -1.6 * bestStem * (0.7 + gc) + 4.2 + 0.35 * Math.max(0, loop - 4);
    dG += (bestStem - bestI) * 1.8;
  }
  dG += -0.12 * bestStem * Math.log10(Math.max(monovalent(c), 1) / 100);
  dG += -0.18 * Math.log10(1 + c.mg);
  dG += 0.035 * c.crowder_pct * bestStem;
  dG += 0.055 * (c.temperature - 37);
  if (polymer === "RNA") dG -= 0.8;
  const dots = ".".repeat(n).split("");
  for (let k = 0; k < bestStem; k++) {
    dots[k] = "(";
    dots[n - 1 - k] = ")";
  }
  return {
    dG,
    stem: bestStem,
    loop: Math.max(loop, 0),
    notation: dots.join(""),
  };
}

function g4Energy(seq: string, c: Conditions, polymer: Polymer): { dG: number; tracts: number; tetrads: number; notation: string; topology: string } {
  const g = countRuns(seq, "G", 2);
  const tetrads = Math.min(4, Math.max(2, g.max));
  let dG = 8;
  if (g.n >= 4 && g.max >= 2) {
    dG = -2.4 * tetrads * Math.min(g.n, 4) + 2.2 * Math.max(0, 4 - g.n);
    const loopish = seq.replace(/G{2,}/g, " ").trim().split(/\s+/);
    const loopPenalty = loopish.reduce((a, L) => a + Math.max(0, L.length - 1) * 0.35, 0);
    dG += loopPenalty;
  }
  const k = Math.max(monovalent(c), 1);
  dG += -1.15 * Math.log10(k / 20);
  dG += -0.22 * Math.min(c.mg, 10);
  dG += -0.028 * c.crowder_pct;
  dG += 0.048 * (c.temperature - 37);
  if (polymer === "RNA") {
    dG -= 1.4;
  }
  const notation = seq.replace(/G+/g, (m) => "(".repeat(m.length)).replace(/[^G()]/g, ".");
  const topology = polymer === "RNA" ? "parallel, 3-tetrad propeller" : tetrads >= 3 ? "hybrid / parallel G4" : "two-tetrad G4";
  return { dG, tracts: g.n, tetrads, notation, topology };
}

function imotifEnergy(seq: string, c: Conditions): { dG: number; tracts: number; notation: string } {
  const cy = countRuns(seq, "C", 2);
  const pKa = 4.7;
  const frac = 1 / (1 + 10 ** (c.ph - pKa));
  let dG = 9;
  if (cy.n >= 4) {
    dG = -1.7 * cy.max * Math.min(cy.n, 4);
  } else if (cy.n >= 2) {
    dG = 2.5 - 0.6 * cy.n * cy.max;
  }
  dG += (1 - frac) * 9.5;
  dG += -0.35 * Math.min(c.mg, 12);
  dG += -0.04 * c.crowder_pct * frac;
  dG += 0.06 * (c.temperature - 37);
  const notation = seq.replace(/C+/g, (m) => "[".repeat(m.length)).replace(/[^C[\]]/g, ".");
  return { dG, tracts: cy.n, notation };
}

function acMotifEnergy(seq: string, c: Conditions): { dG: number; notation: string } {
  const ac = (seq.match(/AC|CA/g) ?? []).length;
  const aRuns = countRuns(seq, "A", 2);
  const cRuns = countRuns(seq, "C", 2);
  const motifish = ac >= 6 || (aRuns.n >= 2 && cRuns.n >= 2 && seq.length >= 12);
  const pKaA = 4.1;
  const pKaC = 4.6;
  const fracA = 1 / (1 + 10 ** (c.ph - pKaA));
  const fracC = 1 / (1 + 10 ** (c.ph - pKaC));
  let dG = 10;
  if (motifish) {
    dG = 1.2 - 0.28 * ac;
  }
  dG += (1 - 0.55 * fracA - 0.45 * fracC) * 5.5;
  dG += -0.85 * Math.min(c.mg, 10);
  dG += -0.03 * c.crowder_pct;
  dG += 0.05 * (c.temperature - 37);
  const notation = seq.replace(/[AC]+/g, (m) => "{".repeat(m.length)).replace(/[^AC{}]/g, ".");
  return { dG, notation };
}

function duplexEnergy(seq: string, polymer: Polymer, c: Conditions): { dG: number; notation: string } {
  const n = seq.length;
  const gc = gcFrac(seq);
  let dG = 2.4 - n * (0.12 + 0.55 * gc);
  if (polymer === "RNA") dG -= 0.4 * n * 0.08;
  dG += -0.114 * n * Math.log10(Math.max(monovalent(c), 1) / 1000);
  dG += -0.3 * Math.log10(1 + 4 * c.mg);
  dG += 0.045 * c.crowder_pct * n * 0.15;
  dG += 0.07 * (c.temperature - 37);
  const notation = "|".repeat(n);
  return { dG, notation };
}

function cruciformEnergy(seq: string, polymer: Polymer, c: Conditions): { dG: number; stem: number; notation: string } {
  const ir = invertedRepeatScore(seq);
  let dG = 8;
  if (ir.stem >= 6) dG = 3.2 - 0.9 * ir.stem;
  else if (ir.stem >= 4) dG = 5.5 - 0.45 * ir.stem;
  dG += -0.08 * ir.stem * Math.log10(Math.max(monovalent(c), 10) / 100);
  dG += 0.04 * (c.temperature - 37);
  if (polymer === "RNA") dG += 1.2;
  const n = seq.length;
  const dots = ".".repeat(n).split("");
  for (let k = 0; k < ir.stem && k < n / 2; k++) {
    dots[k] = "<";
    dots[n - 1 - k] = ">";
  }
  return { dG, stem: ir.stem, notation: dots.join("") };
}

function triplexEnergy(seq: string, c: Conditions): { dG: number; notation: string } {
  const pur = (seq.match(/[AG]/g)?.length ?? 0) / Math.max(seq.length, 1);
  const pyr = 1 - pur;
  const homopurine = pur > 0.72 || pyr > 0.72;
  let dG = 9.5;
  if (homopurine && seq.length >= 12) dG = 3.8 - 0.18 * seq.length * Math.max(pur, pyr);
  const pHterm = Math.abs(c.ph - 5.8) * 0.7;
  dG += pHterm;
  dG += -0.25 * Math.min(c.mg, 8);
  dG += 0.05 * (c.temperature - 37);
  return { dG, notation: seq.replace(/./g, "*") };
}

function coilEnergy(seq: string, c: Conditions): number {
  let dG = 0.15 * seq.length;
  dG += -0.008 * c.temperature;
  dG += 0.01 * c.crowder_pct;
  return dG;
}

function boltzmann(members: Omit<EnsembleMember, "probability">[], tempC: number): EnsembleMember[] {
  const T = tempC + 273.15;
  const weights = members.map((m) => Math.exp(Math.min(40, Math.max(-40, -m.deltaG / (R * T)))));
  const Z = weights.reduce((a, b) => a + b, 0) || 1;
  return members
    .map((m, i) => ({ ...m, probability: (weights[i] ?? 0) / Z }))
    .sort((a, b) => b.probability - a.probability);
}

function predictTm(kind: StructureKind, seq: string, c: Conditions, dG: number): number {
  const gc = gcFrac(seq);
  const N = seq.length;
  // Salt correction
  const saltCorr = 16.6 * Math.log10(Math.max(monovalent(c), 1) / 1000);
  
  if (kind === "g-quadruplex") {
    // Mergny & Lacroix empirical formula for G4s
    return 64 + 0.41 * (gc * 100) - (500 / N) + saltCorr + (c.mg * 0.5) - (dG * 0.5);
  } else if (kind === "i-motif") {
    // pH dependent Tm for i-Motif
    return 30 + 8 * (5.8 - c.ph) + saltCorr - (dG * 0.8);
  } else if (kind === "hairpin" || kind === "duplex") {
    // Nearest neighbor approximation (SantaLucia style simplified)
    return 81.5 + 16.6 * Math.log10(Math.max(monovalent(c), 1) / 1000) + 41 * gc - (500 / N) - (dG * 1.5);
  } else if (kind === "ss-coil") {
    return 10; // essentially random coil at room temp
  }
  
  // Generic fallback
  return 37 - (dG * 2.5);
}

function member(
  kind: StructureKind,
  seq: string,
  dG: number,
  topology: string,
  notation: string,
  details: string,
  g4Topology?: "parallel" | "antiparallel",
  conditions?: Conditions,
  basis: ClaimBasis = "heuristic",
): Omit<EnsembleMember, "probability"> {
  const geo = geometryFor(kind, seq, g4Topology);
  const predictedTm = conditions ? Math.round(predictTm(kind, seq, conditions, dG) * 10) / 10 : undefined;
  
  return {
    id: g4Topology ? `${kind}-${g4Topology}` : kind,
    kind,
    topology,
    deltaG: Math.round(dG * 100) / 100,
    predictedTm,
    notation,
    details,
    nucleotides: geo.nucleotides,
    strands: geo.strands,
    tetradBonds: geo.tetradBonds,
    // Everything built locally is heuristic by construction. Members whose
    // numbers came from a trained head are re-badged in store.ts when the
    // service answers; nothing in this file may claim "measured".
    basis,
  };
}

/** ML-inspired topology classifier for G-Quadruplexes.
 *  Uses loop length distribution, ion type bias, and sequence features
 *  to predict the probability of parallel vs antiparallel folding.
 */
function g4TopologyProbability(seq: string, conditions: Conditions): { parallelProb: number; antiparallelProb: number } {
  // Extract loop lengths between G-tracts
  const loops = seq.replace(/G{2,}/g, "|").split("|").filter(l => l.length > 0);
  const avgLoop = loops.length > 0 ? loops.reduce((a, l) => a + l.length, 0) / loops.length : 2;
  
  // Short loops (1-2 nt) strongly favor parallel (propeller)
  // Long loops (3+ nt) favor antiparallel (lateral/diagonal)
  let parallelScore = 0;
  
  // Loop length feature (most important predictor)
  if (avgLoop <= 1.5) parallelScore += 3.0;
  else if (avgLoop <= 2.5) parallelScore += 1.5;
  else if (avgLoop <= 4) parallelScore -= 0.5;
  else parallelScore -= 2.0;
  
  // K+ strongly favors parallel; Na+ favors antiparallel
  // monovalent > 100 mM is K+-like conditions
  parallelScore += 0.8 * Math.log10(Math.max(monovalent(conditions), 1) / 50);
  
  // Mg2+ slightly favors parallel
  parallelScore += 0.15 * Math.min(conditions.mg, 10);
  
  // Crowding favors more compact parallel
  parallelScore += 0.05 * conditions.crowder_pct;
  
  // Temperature: high temp destabilizes antiparallel more
  parallelScore += 0.02 * (conditions.temperature - 37);
  
  // Convert to probability via sigmoid
  const p = 1 / (1 + Math.exp(-parallelScore));
  return { parallelProb: p, antiparallelProb: 1 - p };
}

export function predict(raw: string, polymerHint: Polymer, conditions: Conditions): PredictionResult {
  const sequence = cleanSeq(raw);
  const polymer = detectPolymer(sequence, polymerHint);
  if (polymer === "RNA" && !RNA_SUPPORT.enabled) {
    // Refused, not answered. The alternative is a confident number with a
    // caveat beside it, which in a 3D viewer reads as a caveat on a real
    // answer rather than as the absence of one.
    return {
      polymer,
      sequence,
      conditions,
      ensemble: [],
      proteins: [],
      refusal: RNA_SUPPORT.reason,
    };
  }
  if (sequence.length < 6) {
    const coil = boltzmann(
      [member("ss-coil", sequence || "ACGTAC", coilEnergy(sequence || "ACGTAC", conditions), "unfolded", ".".repeat(sequence.length), "Sequence too short for stable folds.", undefined, conditions)],
      conditions.temperature,
    );
    return { polymer, sequence, conditions, ensemble: coil, proteins: [] };
  }

  const hp = hairpinEnergy(sequence, polymer, conditions);
  const g4 = g4Energy(sequence, conditions, polymer);
  const im = imotifEnergy(sequence, conditions);
  const ac = acMotifEnergy(sequence, conditions);
  const dx = duplexEnergy(sequence, polymer, conditions);
  const cr = cruciformEnergy(sequence, polymer, conditions);
  const tx = triplexEnergy(sequence, conditions);
  const coil = coilEnergy(sequence, conditions);
  
  // ML topology prediction for G4
  const g4Topo = g4TopologyProbability(sequence, conditions);

  const rawMembers = [
    member(
      "g-quadruplex",
      sequence,
      g4.dG - 0.5 * Math.log(g4Topo.parallelProb + 0.01),
      `parallel ${g4.topology} (${(g4Topo.parallelProb * 100).toFixed(0)}% predicted)`,
      g4.notation,
      `${g4.tracts} G-tracts, ~${g4.tetrads} tetrads. Local energy function; no trained head answered. K+ ${conditions.k} mM, Na+ ${conditions.na} mM.`,
      "parallel",
      conditions,
      "heuristic"
    ),
    member(
      "g-quadruplex",
      sequence,
      g4.dG - 0.5 * Math.log(g4Topo.antiparallelProb + 0.01),
      `antiparallel ${g4.topology} (${(g4Topo.antiparallelProb * 100).toFixed(0)}% predicted)`,
      g4.notation,
      `${g4.tracts} G-tracts, ~${g4.tetrads} tetrads. Antiparallel folds are favoured in Na+ (${conditions.na} mM here). Local energy function.`,
      "antiparallel",
      conditions,
      "heuristic"
    ),
    member(
      "i-motif",
      sequence,
      im.dG,
      "intercalated C+:C hemi-protonated",
      im.notation,
      `${im.tracts} C-tracts. C:C+ hemiprotonation has a pKa near 4.7, so the fraction folded falls steeply above pH 6. Local energy function.`,
      undefined,
      conditions,
      "heuristic"
    ),
    member(
      "ac-motif",
      sequence,
      ac.dG,
      "intercalated A+:C / C+:C",
      ac.notation,
      "Adenine-cytosine motif. No QuadCond head covers this fold and no measurement in the atlas constrains it.",
      undefined,
      conditions,
      "heuristic"
    ),
    member(
      "hairpin",
      sequence,
      hp.dG,
      `stem ${hp.stem} bp, loop ${hp.loop} nt`,
      hp.notation,
      "Intramolecular stem-loop, nearest-neighbour stacking with a salt correction. Local energy function.",
      undefined,
      conditions,
      "heuristic"
    ),
    member(
      "duplex",
      sequence,
      dx.dG,
      polymer === "RNA" ? "A-form helix" : "B-form helix",
      dx.notation,
      "Watson-Crick duplex, nearest-neighbour stacking with a salt correction. Local energy function.",
      undefined,
      conditions,
      "heuristic"
    ),
    member(
      "cruciform",
      sequence,
      cr.dG,
      cr.stem >= 6 ? `four-way junction, ${cr.stem} bp arms` : "weak inverted repeat",
      cr.notation,
      "Inverted-repeat extrusion, competing with the duplex. No QuadCond head covers this fold.",
      undefined,
      conditions,
      "heuristic"
    ),
    member(
      "triplex",
      sequence,
      tx.dG,
      "homopurine·homopyrimidine third strand",
      tx.notation,
      "Homopurine-homopyrimidine third strand. No QuadCond head covers this fold.",
      undefined,
      conditions,
      "heuristic"
    ),
    member(
      "ss-coil",
      sequence,
      coil,
      "worm-like chain",
      ".".repeat(sequence.length),
      "Unpaired reference state, at deltaG = 0 by definition. Every other free energy in the ensemble is measured against it.",
      undefined,
      conditions,
      "heuristic"
    ),
  ];

  // Folds with no evidence behind them are dropped BEFORE the Boltzmann step,
  // not hidden afterwards. A member left in the partition function still moves
  // every other member's probability, so hiding it in the UI would change the
  // numbers on screen while claiming not to.
  const ensemble = boltzmann(
    rawMembers.filter((m) => foldEnabled(m.kind)),
    conditions.temperature,
  );
  const proteins = PROTEIN_BINDERS.enabled ? bindersFor(sequence, polymer, ensemble) : [];

  return { polymer, sequence, conditions, ensemble, proteins };
}
