import { type EnsembleMember } from "./types";
import { MIN_CONTACT, MODEL_DISCLAIMER } from "./geometry-params";

/**
 * Write a coarse-grained PDB model of an ensemble member.
 *
 * Three things were wrong with the previous version and each of them produces a
 * file that a structure viewer either rejects or silently misreads:
 *
 *  - every atom was named ` P  ` and given element `P`, so a 22-nucleotide DNA
 *    fold was written as 22 phosphorus atoms and rendered as a cloud of
 *    unconnected dots;
 *  - residue names were the raw one-letter bases (`A`, `G`), which are amino
 *    acid codes in PDB — a viewer reads `G` as glycine;
 *  - CONECT serials were derived from `nt.index`, the position in the
 *    *sequence*, while ATOM serials came from a running counter. Any structure
 *    whose residues were not emitted in sequence order got bonds between the
 *    wrong atoms.
 *
 * What is written now is an explicit coarse-grained model: one C1' pseudo-atom
 * per residue, standard DNA/RNA residue names, one chain per strand, and a
 * REMARK block that says in the file itself that these are idealised model
 * coordinates. That last part matters more than the rest — a file that does not
 * announce its own status will be read as a prediction of atomic coordinates,
 * which it is not.
 */

const DNA_RESIDUE: Record<string, string> = {
  A: " DA", C: " DC", G: " DG", T: " DT", U: " DU", N: " DN",
};
const RNA_RESIDUE: Record<string, string> = {
  A: "  A", C: "  C", G: "  G", U: "  U", T: "  T", N: "  N",
};

const CHAIN_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

function fixed(value: number, width: number, decimals: number): string {
  return value.toFixed(decimals).padStart(width, " ");
}

export interface PdbOptions {
  polymer?: "DNA" | "RNA";
  sequenceId?: string;
  conditions?: string;
  probability?: number;
}

export function buildPDB(member: EnsembleMember, opts: PdbOptions = {}): string {
  const polymer = opts.polymer ?? "DNA";
  const table = polymer === "RNA" ? RNA_RESIDUE : DNA_RESIDUE;
  const lines: string[] = [];

  lines.push(`HEADER    ${member.kind.toUpperCase().padEnd(40)}`);
  lines.push(`TITLE     AENNA-3D MODEL OF ${(opts.sequenceId ?? member.id).toUpperCase()}`);
  lines.push("REMARK   1");
  lines.push("REMARK   1 THIS IS A MODEL, NOT AN EXPERIMENTAL STRUCTURE.");
  for (const chunk of MODEL_DISCLAIMER.match(/.{1,66}(\s|$)/g) ?? []) {
    lines.push(`REMARK   1 ${chunk.trim().toUpperCase()}`);
  }
  lines.push("REMARK   1");
  lines.push("REMARK   2 RESOLUTION. NOT APPLICABLE.");
  lines.push("REMARK 220 COARSE-GRAINED: ONE C1' PSEUDO-ATOM PER RESIDUE.");
  lines.push("REMARK 220 SIDE CHAIN, SUGAR PUCKER AND BACKBONE TORSIONS ARE NOT MODELLED.");
  lines.push(`REMARK 220 TOPOLOGY: ${(member.topology ?? "N/A").toUpperCase()}`);
  if (opts.conditions) lines.push(`REMARK 220 CONDITIONS: ${opts.conditions.toUpperCase()}`);
  if (typeof member.deltaG === "number") {
    lines.push(`REMARK 220 MODEL DELTA-G (KCAL/MOL): ${member.deltaG.toFixed(2)}`);
  }
  if (typeof opts.probability === "number") {
    lines.push(`REMARK 220 ENSEMBLE WEIGHT: ${opts.probability.toFixed(4)}`);
  }
  lines.push("MODEL        1");

  // ATOM serial must map 1:1 to the array position, because CONECT refers to
  // serials and the strand arrays refer to array positions.
  const serialOf = new Map<number, number>();
  const nts = member.nucleotides;
  let serial = 1;
  for (let i = 0; i < nts.length; i++) {
    const nt = nts[i];
    if (!nt) continue;
    const resName = table[nt.base.toUpperCase()] ?? table.N!;
    const chainId = CHAIN_IDS[nt.strand % CHAIN_IDS.length] ?? "A";
    const bFactor = nt.role === "tetrad" || nt.role === "pair" ? 20.0 : 60.0;
    lines.push(
      "ATOM  " +
        serial.toString().padStart(5, " ") +
        " " + " C1'" +
        " " + resName +
        " " + chainId +
        (nt.index + 1).toString().padStart(4, " ") +
        "    " +
        fixed(nt.x, 8, 3) + fixed(nt.y, 8, 3) + fixed(nt.z, 8, 3) +
        fixed(1.0, 6, 2) + fixed(bFactor, 6, 2) +
        "          " + " C",
    );
    serialOf.set(i, serial);
    serial++;
  }
  lines.push("TER   " + serial.toString().padStart(5, " "));
  lines.push("ENDMDL");

  // Backbone connectivity, then the H-bond network that defines the fold.
  for (const strandIdx of member.strands) {
    for (let i = 0; i < strandIdx.length - 1; i++) {
      const a = serialOf.get(strandIdx[i] as number);
      const b = serialOf.get(strandIdx[i + 1] as number);
      if (a && b) {
        lines.push("CONECT" + a.toString().padStart(5, " ") + b.toString().padStart(5, " "));
      }
    }
  }
  for (const [i, j] of member.tetradBonds ?? []) {
    const a = serialOf.get(i);
    const b = serialOf.get(j);
    if (a && b) {
      lines.push("CONECT" + a.toString().padStart(5, " ") + b.toString().padStart(5, " "));
    }
  }
  lines.push("END");
  return lines.join("\n") + "\n";
}

/** Geometry self-check. Returns the problems found, empty when the model is sane. */
export function checkGeometry(member: EnsembleMember): string[] {
  const problems: string[] = [];
  const nts = member.nucleotides;
  const core = nts.filter((n) => n.role === "tetrad" || n.role === "pair");
  for (let i = 0; i < core.length; i++) {
    for (let j = i + 1; j < core.length; j++) {
      const a = core[i]!, b = core[j]!;
      const d = Math.hypot(a.x - b.x, a.y - b.y, a.z - b.z);
      if (d < MIN_CONTACT) {
        problems.push(
          `residues ${a.index + 1} and ${b.index + 1} are ${d.toFixed(2)} A apart ` +
            `(hard-sphere floor ${MIN_CONTACT} A)`,
        );
      }
    }
  }
  if (!nts.every((n) => Number.isFinite(n.x) && Number.isFinite(n.y) && Number.isFinite(n.z))) {
    problems.push("non-finite coordinates");
  }
  return problems;
}

export function exportToPDB(member: EnsembleMember, filename = "structure.pdb",
                            opts: PdbOptions = {}) {
  const text = buildPDB(member, opts);
  const blob = new Blob([text], { type: "chemical/x-pdb" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}
