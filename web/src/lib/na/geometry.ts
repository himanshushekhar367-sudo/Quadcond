import type { Nucleotide, StructureKind } from "./types";
import { BACKBONE_STEP, DUPLEX, G4, IMOTIF } from "./geometry-params";

// B-form duplex. These were the only correct constants in the original file,
// and the bug was that every other fold was built from them too.
const RISE = DUPLEX.rise;
const TWIST = (DUPLEX.twistDeg * Math.PI) / 180;
const HELIX_R = DUPLEX.c1Radius;

function clean(seq: string): string {
  return seq.toUpperCase().replace(/[^ACGTU]/g, "");
}

function helixNt(
  seq: string,
  i: number,
  strand: number,
  y0: number,
  phase: number,
): Nucleotide {
  const ang = i * TWIST + phase + (strand === 1 ? Math.PI * 0.86 : 0);
  const r = HELIX_R;
  return {
    index: i,
    base: seq[i] ?? "N",
    x: r * Math.cos(ang),
    y: y0 + i * RISE,
    z: r * Math.sin(ang),
    paired: true,
    role: "stem",
    strand,
  };
}

export function buildDuplex(seq: string): { nucleotides: Nucleotide[]; strands: number[][] } {
  const s = clean(seq);
  const nts: Nucleotide[] = [];
  const st0: number[] = [];
  const st1: number[] = [];
  for (let i = 0; i < s.length; i++) {
    nts.push(helixNt(s, i, 0, 0, 0));
    st0.push(nts.length - 1);
  }
  const rc = reverseComplement(s);
  for (let i = 0; i < rc.length; i++) {
    const nt = helixNt(rc, s.length - 1 - i, 1, 0, 0);
    nt.index = i;
    nt.base = rc[i] ?? "N";
    nts.push(nt);
    st1.push(nts.length - 1);
  }
  return { nucleotides: nts, strands: [st0, st1] };
}

export function buildHairpin(seq: string): { nucleotides: Nucleotide[]; strands: number[][] } {
  const s = clean(seq);
  const n = s.length;
  const stem = Math.max(3, Math.floor((n - 4) / 2));
  const loopN = Math.max(3, n - 2 * stem);
  const nts: Nucleotide[] = [];
  const strand: number[] = [];

  for (let i = 0; i < stem; i++) {
    const nt = helixNt(s, i, 0, 0, 0);
    nt.role = "stem";
    nts.push(nt);
    strand.push(nts.length - 1);
  }
  const topY = stem * RISE + 4;
  for (let k = 0; k < loopN; k++) {
    const t = (k + 1) / (loopN + 1);
    const ang = Math.PI * t;
    nts.push({
      index: stem + k,
      base: s[stem + k] ?? "N",
      x: HELIX_R * Math.cos(ang * 0.35),
      y: topY + Math.sin(ang) * 7.5,
      z: HELIX_R * Math.sin(Math.PI * 0.15) + Math.cos(ang) * 6,
      paired: false,
      role: "loop",
      strand: 0,
    });
    strand.push(nts.length - 1);
  }
  for (let i = 0; i < stem; i++) {
    const idx = stem + loopN + i;
    const back = stem - 1 - i;
    const nt = helixNt(s, back, 1, 0, 0);
    nt.index = idx;
    nt.base = s[idx] ?? "N";
    nt.role = "stem";
    nts.push(nt);
    strand.push(nts.length - 1);
  }
  return { nucleotides: nts, strands: [strand] };
}

/**
 * The three G-quadruplex topologies, and which way each column runs.
 *
 * - **parallel** — all four tracts in the same direction (propeller loops).
 * - **antiparallel** — adjacent columns alternate (lateral / diagonal loops).
 * - **hybrid (3+1)** — three columns one way, the fourth the other. This is the
 *   form human telomeric DNA adopts in K⁺, so it is not an exotic case: it is
 *   what the shipped hTelo examples actually predict. It was missing entirely,
 *   and a hybrid card resolved to the parallel archetype — a different fold,
 *   drawn confidently, with nothing on screen saying so.
 */
export type G4Topology = "parallel" | "antiparallel" | "hybrid";

function columnGoesUp(topology: G4Topology, column: number): boolean {
  if (topology === "parallel") return true;
  if (topology === "hybrid") return column !== 3; // 3 + 1
  return column % 2 === 0;
}

export function buildG4(seq: string, topology: G4Topology = "parallel"): { nucleotides: Nucleotide[]; strands: number[][]; tetradBonds?: [number, number][] } {
  const s = clean(seq);
  const nts: Nucleotide[] = [];
  const strand: number[] = [];
  const tetradBonds: [number, number][] = [];
  
  // Quadruplex geometry, not duplex geometry. The previous values (r = 8.0,
  // rise = 5.0) were display constants with no relation to a G4: the real
  // stacking distance between tetrads is 3.3 A and the C1' atoms sit 9.7 A from
  // the four-fold axis. See geometry-params.ts.
  const r = G4.c1Radius;
  const rise = G4.rise;
  const twistPerTetrad = (G4.twistDeg * Math.PI) / 180;

  // ---- Parse actual G-tracts from the sequence ----
  const tractRe = /G{2,}/g;
  const tracts: { start: number; end: number; len: number }[] = [];
  let match: RegExpExecArray | null;
  while ((match = tractRe.exec(s)) !== null) {
    tracts.push({ start: match.index, end: match.index + match[0].length, len: match[0].length });
  }

  const nTracts = Math.min(4, tracts.length);
  if (nTracts < 4) {
    // Not enough G-tracts — fall back to coil
    let x = 0, y = 0, z = 0, heading = 0.4;
    for (let i = 0; i < s.length; i++) {
      heading += Math.sin(i * 0.7) * 0.45 + 0.18;
      x += 3.4 * Math.cos(heading);
      z += 3.4 * Math.sin(heading);
      y += 2;
      nts.push({ index: i, base: s[i] ?? "N", x, y, z, paired: false, role: "ss", strand: 0 });
      strand.push(i);
    }
    return { nucleotides: nts, strands: [strand] };
  }

  const t0 = tracts[0]!, t1 = tracts[1]!, t2 = tracts[2]!, t3 = tracts[3]!;
  const levels = Math.min(t0.len, t1.len, t2.len, t3.len);
  
  // Extract loop sequences between tracts
  const loop1Seq = s.slice(t0.end, t1.start);
  const loop2Seq = s.slice(t1.end, t2.start);
  const loop3Seq = s.slice(t2.end, t3.start);
  const headSeq = s.slice(0, t0.start);
  const tailSeq = s.slice(t3.start + levels);

  // Pre-compute grid positions
  const grid: { x: number; y: number; z: number }[][] = [];
  for (let c = 0; c < 4; c++) {
    grid[c] = [];
    for (let lev = 0; lev < levels; lev++) {
      // Successive tetrads are rotated, not stacked in register. Without this
      // term the four columns ran as straight vertical lines and the model had
      // no helical character at all — a quadruplex is right-handed, ~30 deg per
      // tetrad.
      const a = c * (Math.PI / 2) + Math.PI / 4 + lev * twistPerTetrad;
      grid[c]![lev] = { x: r * Math.cos(a), y: lev * rise, z: r * Math.sin(a) };
    }
  }

  const gridIdx: number[][] = Array.from({ length: 4 }, () => Array(levels).fill(-1));
  let seqI = 0;

  // Head residues
  for (let i = 0; i < headSeq.length; i++) {
    nts.push({ index: seqI, base: headSeq[i] ?? "N", x: r + 5 + i * 3, y: 0, z: 0, paired: false, role: "ss", strand: 0 });
    strand.push(nts.length - 1);
    seqI++;
  }

  // Place 4 columns with actual loops between them
  const tractInfos = [t0, t1, t2, t3];
  const loopSeqs = [loop1Seq, loop2Seq, loop3Seq];

  for (let c = 0; c < 4; c++) {
    const goingUp = columnGoesUp(topology, c);

    // Place column nucleotides
    for (let lev = 0; lev < levels; lev++) {
      const level = goingUp ? lev : (levels - 1 - lev);
      const pos = grid[c]![level]!;

      nts.push({
        index: seqI, base: s[seqI] ?? "G",
        x: pos.x, y: pos.y, z: pos.z,
        paired: true, role: "tetrad", strand: 0,
      });
      strand.push(nts.length - 1);
      gridIdx[c]![level] = nts.length - 1;
      seqI++;
    }
    // Advance seqI past any extra G's in the tract
    seqI = tractInfos[c]!.end;

    // Place loop to next column
    if (c < 3) {
      const loopBases = loopSeqs[c]!;
      const nextC = c + 1;
      const nextGoingUp = columnGoesUp(topology, nextC);
      const fromPos = grid[c]![goingUp ? levels - 1 : 0]!;
      const toPos = grid[nextC]![nextGoingUp ? 0 : levels - 1]!;

      for (let k = 0; k < loopBases.length; k++) {
        const t = (k + 1) / (loopBases.length + 1);
        const bulge = 6 * Math.sin(Math.PI * t);
        const midX = (fromPos.x + toPos.x) / 2;
        const midZ = (fromPos.z + toPos.z) / 2;
        const dist = Math.sqrt(midX * midX + midZ * midZ);
        const outX = dist > 0.1 ? midX / dist : 1;
        const outZ = dist > 0.1 ? midZ / dist : 0;

        nts.push({
          index: seqI, base: loopBases[k] ?? "N",
          x: fromPos.x * (1 - t) + toPos.x * t + outX * bulge,
          y: fromPos.y * (1 - t) + toPos.y * t + bulge * 0.5,
          z: fromPos.z * (1 - t) + toPos.z * t + outZ * bulge,
          paired: false, role: "loop", strand: 0,
        });
        strand.push(nts.length - 1);
        seqI++;
      }
    }
  }

  // Tetrad bonds (horizontal Hoogsteen H-bonds)
  for (let lev = 0; lev < levels; lev++) {
    for (let c = 0; c < 4; c++) {
      const a = gridIdx[c]![lev]!;
      const b = gridIdx[(c + 1) % 4]![lev]!;
      if (a >= 0 && b >= 0) tetradBonds.push([a, b]);
    }
  }

  // Tail residues
  seqI = t3.start + levels;
  while (seqI < s.length) {
    const prev = nts[nts.length - 1];
    nts.push({
      index: seqI, base: s[seqI] ?? "N",
      x: (prev?.x ?? 0) + 3, y: (prev?.y ?? 0) - 2, z: prev?.z ?? 0,
      paired: false, role: "ss", strand: 0,
    });
    strand.push(nts.length - 1);
    seqI++;
  }

  return { nucleotides: nts, strands: [strand], tetradBonds };
}

export function buildImotif(
  seq: string,
  kind: "i-motif" | "ac-motif",
): { nucleotides: Nucleotide[]; strands: number[][]; tetradBonds?: [number, number][] } {
  const s = clean(seq);
  const nts: Nucleotide[] = [];
  const strand: number[] = [];
  const hBonds: [number, number][] = [];
  
  // Intercalated geometry. The two rises are not interchangeable: successive
  // C:C+ pairs of ONE duplex are 6.2 A apart because the partner duplex's pairs
  // sit between them, while the intercalated stack advances 3.1 A per step. The
  // previous single 4.5 A rise had nowhere to put the intercalated pairs and
  // made the fold about twice too tall. See geometry-params.ts.
  const r = IMOTIF.c1Radius;
  const rise = IMOTIF.pairRise;
  const duplexRise = IMOTIF.duplexRise;
  const shift = IMOTIF.duplexRise / 2; // half a step: this IS the intercalation
  // `rise` is the height the intercalated stack gains per C:C+ pair counting
  // both duplexes; `duplexRise` is the spacing within one of them. shift ===
  // rise by construction, and naming both makes the relationship checkable
  // rather than a coincidence of two literals.
  void rise;

  // ---- Parse the actual C-tracts (or AC-tracts) and loops from the sequence ----
  const tractBase = kind === "ac-motif" ? "[AC]" : "C";
  const tractRe = new RegExp(`${tractBase}{2,}`, "g");
  const tracts: { start: number; end: number; len: number }[] = [];
  let m: RegExpExecArray | null;
  while ((m = tractRe.exec(s)) !== null) {
    tracts.push({ start: m.index, end: m.index + m[0].length, len: m[0].length });
  }
  
  // We need exactly 4 tracts for an i-Motif; take the first 4
  const nTracts = Math.min(4, tracts.length);
  if (nTracts < 4) {
    // Not enough tracts — fall back to coil-like placement
    let x = 0, y = 0, z = 0, heading = 0.4;
    for (let i = 0; i < s.length; i++) {
      heading += Math.sin(i * 0.7) * 0.45 + 0.18;
      x += 3.4 * Math.cos(heading);
      z += 3.4 * Math.sin(heading);
      y += 2;
      nts.push({ index: i, base: s[i] ?? "N", x, y, z, paired: false, role: "ss", strand: 0 });
      strand.push(i);
    }
    return { nucleotides: nts, strands: [strand] };
  }

  const t0 = tracts[0]!, t1 = tracts[1]!, t2 = tracts[2]!, t3 = tracts[3]!;
  // Use the minimum tract length so all columns are the same height
  const levels = Math.min(t0.len, t1.len, t2.len, t3.len);

  // Extract loop sequences between tracts
  const loop1Seq = s.slice(t0.end, t1.start);
  const loop2Seq = s.slice(t1.end, t2.start);
  const loop3Seq = s.slice(t2.end, t3.start);
  // Head/tail residues outside the 4 tracts
  const headSeq = s.slice(0, t0.start);
  const tailSeq = s.slice(t3.start + levels);

  // Track indices for H-bond pairing
  const colIndices: number[][] = [[], [], [], []];
  let seqI = 0;

  // Place head residues (before first tract)
  for (let i = 0; i < headSeq.length; i++) {
    nts.push({ index: seqI, base: headSeq[i] ?? "N", x: r + 4 + i * 3, y: 0, z: 0, paired: false, role: "ss", strand: 0 });
    strand.push(nts.length - 1);
    seqI++;
  }

  // Col 0: Tract 1 going UP at +X
  for (let p = 0; p < levels; p++) {
    nts.push({ index: seqI, base: s[seqI] ?? "C", x: r, y: p * duplexRise, z: 0, paired: true, role: "pair", strand: 0 });
    strand.push(nts.length - 1);
    colIndices[0]!.push(nts.length - 1);
    seqI++;
  }
  // Skip extra tract bases beyond levels
  seqI = t0.end;

  // Loop 1: top of Col 0 → top of Col 1
  const topY_A = (levels - 1) * duplexRise;
  const topY_B = (levels - 1) * duplexRise + shift;
  for (let k = 0; k < loop1Seq.length; k++) {
    const fromNt = nts[nts.length - 1]!;
    const t = (k + 1) / (loop1Seq.length + 1);
    const bulge = 5 * Math.sin(Math.PI * t);
    nts.push({
      index: seqI, base: loop1Seq[k] ?? "N",
      x: fromNt.x * (1 - t) + 0 * t + bulge * 0.3,
      y: topY_A * (1 - t) + topY_B * t + bulge * 0.4,
      z: fromNt.z * (1 - t) + r * t + bulge * 0.3,
      paired: false, role: "loop", strand: 0,
    });
    strand.push(nts.length - 1);
    seqI++;
  }

  // Col 1: Tract 2 going DOWN at +Z
  for (let p = levels - 1; p >= 0; p--) {
    nts.push({ index: seqI, base: s[seqI] ?? "C", x: 0, y: p * duplexRise + shift, z: r, paired: true, role: "pair", strand: 0 });
    strand.push(nts.length - 1);
    colIndices[1]!.push(nts.length - 1);
    seqI++;
  }
  seqI = t1.end;

  // Loop 2: bottom of Col 1 → bottom of Col 2
  for (let k = 0; k < loop2Seq.length; k++) {
    const fromNt = nts[nts.length - 1]!;
    const t = (k + 1) / (loop2Seq.length + 1);
    const bulge = 5 * Math.sin(Math.PI * t);
    nts.push({
      index: seqI, base: loop2Seq[k] ?? "N",
      x: fromNt.x * (1 - t) + 0 * t - bulge * 0.3,
      y: shift * (1 - t) + shift * t - bulge * 0.4,
      z: fromNt.z * (1 - t) + (-r) * t - bulge * 0.3,
      paired: false, role: "loop", strand: 0,
    });
    strand.push(nts.length - 1);
    seqI++;
  }

  // Col 2: Tract 3 going UP at -Z (pairs with Col 1)
  for (let p = 0; p < levels; p++) {
    nts.push({ index: seqI, base: s[seqI] ?? "C", x: 0, y: p * duplexRise + shift, z: -r, paired: true, role: "pair", strand: 0 });
    strand.push(nts.length - 1);
    colIndices[2]!.push(nts.length - 1);
    seqI++;
  }
  seqI = t2.end;

  // Loop 3: top of Col 2 → top of Col 3
  for (let k = 0; k < loop3Seq.length; k++) {
    const fromNt = nts[nts.length - 1]!;
    const t = (k + 1) / (loop3Seq.length + 1);
    const bulge = 5 * Math.sin(Math.PI * t);
    nts.push({
      index: seqI, base: loop3Seq[k] ?? "N",
      x: fromNt.x * (1 - t) + (-r) * t - bulge * 0.3,
      y: topY_B * (1 - t) + topY_A * t + bulge * 0.4,
      z: fromNt.z * (1 - t) + 0 * t + bulge * 0.3,
      paired: false, role: "loop", strand: 0,
    });
    strand.push(nts.length - 1);
    seqI++;
  }

  // Col 3: Tract 4 going DOWN at -X (pairs with Col 0)
  for (let p = levels - 1; p >= 0; p--) {
    nts.push({ index: seqI, base: s[seqI] ?? "C", x: -r, y: p * duplexRise, z: 0, paired: true, role: "pair", strand: 0 });
    strand.push(nts.length - 1);
    colIndices[3]!.push(nts.length - 1);
    seqI++;
  }
  seqI = t3.start + levels;

  // Build C·C+ hydrogen bonds
  // Col 0 (up: levels 0,1,2) ↔ Col 3 (down: levels 2,1,0)
  const col0 = colIndices[0]!, col3 = colIndices[3]!;
  for (let p = 0; p < Math.min(col0.length, col3.length); p++) {
    hBonds.push([col0[p]!, col3[col3.length - 1 - p]!]);
  }
  // Col 1 (down: levels 2,1,0) ↔ Col 2 (up: levels 0,1,2)
  const col1 = colIndices[1]!, col2 = colIndices[2]!;
  for (let p = 0; p < Math.min(col1.length, col2.length); p++) {
    hBonds.push([col1[p]!, col2[col2.length - 1 - p]!]);
  }

  // Tail residues
  while (seqI < s.length) {
    const prev = nts[nts.length - 1];
    nts.push({
      index: seqI, base: s[seqI] ?? "N",
      x: (prev?.x ?? 0) - 3, y: (prev?.y ?? 0) - 2, z: (prev?.z ?? 0),
      paired: false, role: "ss", strand: 0,
    });
    strand.push(nts.length - 1);
    seqI++;
  }

  return { nucleotides: nts, strands: [strand], tetradBonds: hBonds };
}

export function buildCruciform(seq: string): { nucleotides: Nucleotide[]; strands: number[][] } {
  const s = clean(seq);
  const nts: Nucleotide[] = [];
  const arms: number[][] = [[], [], [], []];
  const armLen = Math.max(4, Math.floor(s.length / 4));
  const dirs: [number, number][] = [
    [1, 0],
    [0, 1],
    [-1, 0],
    [0, -1],
  ];
  let idx = 0;
  for (let a = 0; a < 4; a++) {
    const [dx, dz] = dirs[a] ?? [1, 0];
    for (let k = 0; k < armLen; k++) {
      const y = (k % 2) * 1.6;
      const r = 4 + k * RISE * 0.55;
      nts.push({
        index: idx,
        base: s[idx] ?? "N",
        x: dx * r,
        y,
        z: dz * r,
        paired: k > 0,
        role: k === 0 ? "junction" : "stem",
        strand: a,
      });
      arms[a]?.push(nts.length - 1);
      idx++;
      if (idx >= s.length) break;
    }
  }
  while (idx < s.length) {
    nts.push({
      index: idx,
      base: s[idx] ?? "N",
      x: 2 + idx * 0.4,
      y: 3,
      z: 2,
      paired: false,
      role: "ss",
      strand: 0,
    });
    arms[0]?.push(nts.length - 1);
    idx++;
  }
  return { nucleotides: nts, strands: arms.filter((x) => x.length) };
}

export function buildTriplex(seq: string): { nucleotides: Nucleotide[]; strands: number[][] } {
  const s = clean(seq);
  const n = Math.ceil(s.length / 3);
  const nts: Nucleotide[] = [];
  const strands: number[][] = [[], [], []];
  const radii = [8.6, 8.6, 11.4];
  const phases = [0, Math.PI * 0.86, Math.PI * 0.43];
  let idx = 0;
  for (let st = 0; st < 3; st++) {
    for (let i = 0; i < n && idx < s.length; i++) {
      const ang = i * TWIST + (phases[st] ?? 0);
      const r = radii[st] ?? 9;
      nts.push({
        index: idx,
        base: s[idx] ?? "N",
        x: r * Math.cos(ang),
        y: i * RISE,
        z: r * Math.sin(ang),
        paired: true,
        role: "stem",
        strand: st,
      });
      strands[st]?.push(nts.length - 1);
      idx++;
    }
  }
  return { nucleotides: nts, strands };
}

export function buildCoil(seq: string): { nucleotides: Nucleotide[]; strands: number[][] } {
  const s = clean(seq);
  const nts: Nucleotide[] = [];
  const strand: number[] = [];
  let x = 0,
    y = 0,
    z = 0;
  let heading = 0.4;
  for (let i = 0; i < s.length; i++) {
    heading += Math.sin(i * 0.7) * 0.45 + 0.18;
    const pitch = Math.sin(i * 0.33) * 0.35;
    x += 3.4 * Math.cos(heading);
    z += 3.4 * Math.sin(heading);
    y += 3.4 * Math.sin(pitch);
    nts.push({
      index: i,
      base: s[i] ?? "N",
      x,
      y,
      z,
      paired: false,
      role: "ss",
      strand: 0,
    });
    strand.push(i);
  }
  return { nucleotides: nts, strands: [strand] };
}

export function geometryFor(
  kind: StructureKind,
  seq: string,
  g4Topology?: G4Topology,
): { nucleotides: Nucleotide[]; strands: number[][]; tetradBonds?: [number, number][] } {
  switch (kind) {
    case "duplex":
      return buildDuplex(seq);
    case "hairpin":
      return buildHairpin(seq);
    case "g-quadruplex":
      return buildG4(seq, g4Topology ?? "parallel");
    case "i-motif":
      return buildImotif(seq, "i-motif");
    case "ac-motif":
      return buildImotif(seq, "ac-motif");
    case "cruciform":
      return buildCruciform(seq);
    case "triplex":
      return buildTriplex(seq);
    default:
      return buildCoil(seq);
  }
}

function reverseComplement(seq: string): string {
  const map: Record<string, string> = {
    A: seq.includes("U") ? "U" : "T",
    T: "A",
    U: "A",
    G: "C",
    C: "G",
  };
  return seq
    .split("")
    .reverse()
    .map((b) => map[b] ?? "N")
    .join("");
}

function findGTracts(seq: string): { nTetrads: number } {
  const m = seq.match(/G{2,6}/g) ?? [];
  const nTetrads = m.reduce((a, t) => a + t.length, 0) >= 12 ? 3 : 2;
  return { nTetrads };
}
