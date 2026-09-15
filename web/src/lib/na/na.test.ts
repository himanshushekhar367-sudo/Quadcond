/**
 * Tests for the parts of the viewer that make claims.
 *
 * Not the rendering — the rendering is judged by looking at it. These cover the
 * places where a wrong answer would look right, and the list is drawn from
 * defects this codebase actually shipped rather than from what seemed worth
 * checking. Two of the earlier tests here were themselves the problem: one
 * "distinguishes potassium from sodium" test asserted only that two input
 * objects differed, and passed happily while the function behind them returned
 * byte-identical output for both.
 */
import { describe, it } from "node:test";
import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";

import { geometryOnly, memberForCard } from "./geometry-only.ts";
import { DEFAULT_CONDITIONS, monovalent, type Conditions } from "./types.ts";
import { FOLD_CAPABILITIES, foldEnabled, disabledFolds } from "./capabilities.ts";
import { visibleFolds, extrapolatedFolds, claimBadge } from "../quadcond/client.ts";
import type { EvidenceCard } from "../quadcond/client.ts";
import { templatesFor, templatesForKind, provenanceLabel, TEMPLATES } from "./templates.ts";

const HTELO = "AGGGTTAGGGTTAGGGTTAGGG";

function conditions(over: Partial<Conditions> = {}): Conditions {
  return { ...DEFAULT_CONDITIONS, ...over };
}

describe("the local path produces no numbers", () => {
  it("returns structures with no deltaG and no probability", () => {
    // The whole point of geometry-only. A local energy function that could not
    // tell K+ from Na+, was never given the dG calibration, and put its own
    // reference state at +3 kcal/mol was rendering probability bars under a
    // banner that said nothing was being substituted.
    const r = geometryOnly(HTELO, "DNA", conditions());
    assert.ok(r.members.length > 0, "structures should still be drawn");
    for (const m of r.members) {
      assert.ok(!("deltaG" in m), `${m.id} carries a deltaG`);
      assert.ok(!("probability" in m), `${m.id} carries a probability`);
      assert.ok(!("predictedTm" in m), `${m.id} carries a Tm`);
    }
    assert.ok(!("ensemble" in r), "the field is gone, not empty — an empty array kept stale reads compiling");
  });

  it("draws the same structures whatever the buffer, because buffer is not its business", () => {
    const a = geometryOnly(HTELO, "DNA", conditions({ k: 150, na: 0, ph: 5 }));
    const b = geometryOnly(HTELO, "DNA", conditions({ k: 0, na: 150, ph: 8 }));
    assert.deepEqual(
      a.members.map((m) => m.id),
      b.members.map((m) => m.id),
    );
    // and crucially there is nothing numeric that could have differed silently
    assert.equal(JSON.stringify(a.members.map((m) => Object.keys(m).sort())),
                 JSON.stringify(b.members.map((m) => Object.keys(m).sort())));
  });

  it("refuses RNA outright rather than drawing it", () => {
    const r = geometryOnly("GGGUUAGGGUUAGGGUUAGGG", "RNA", conditions());
    assert.equal(r.members.length, 0);
    assert.match(r.refusal ?? "", /DNA/);
  });

  it("never draws a fold with no evidence behind it", () => {
    const r = geometryOnly(HTELO, "DNA", conditions());
    for (const kind of ["ac-motif", "cruciform", "triplex"] as const) {
      assert.equal(foldEnabled(kind), false);
      assert.equal(r.members.some((m) => m.kind === kind), false);
    }
  });
});

describe("service cards map onto drawable geometry, by strand", () => {
  // The defect: `geometryOnly` built members only for the typed strand, so the
  // reverse-strand card had no geometry and ViewerPane fell back to
  // `members[0]`. On a G-rich input the panel said "i-motif, strand −" and the
  // picture showed a parallel G4. Two molecules, one screen, no error.
  const G = "AGGGTTAGGGTTAGGGTTAGGG";
  const C = "CCCTAACCCTAACCCTAACCCT";

  it("builds the complementary-strand structure for a G-rich input", () => {
    const r = geometryOnly(G, "DNA", conditions());
    const im = r.members.find((m) => m.kind === "i-motif");
    assert.ok(im, "a G-rich input must still offer the complementary i-motif");
    assert.equal(im.strand, "-");
    assert.ok(im.strandSequence.startsWith("CCC"), im.strandSequence);
  });

  it("builds the complementary-strand structure for a C-rich input", () => {
    const r = geometryOnly(C, "DNA", conditions());
    const g4 = r.members.find((m) => m.kind === "g-quadruplex");
    assert.ok(g4, "a C-rich input must still offer the complementary G4");
    assert.equal(g4.strand, "-");
    assert.ok(g4.strandSequence.startsWith("AGGG"), g4.strandSequence);
  });

  it("matches a card to the member on the SAME strand", () => {
    const r = geometryOnly(G, "DNA", conditions());
    const hit = memberForCard(r.members, { kind: "i-motif", strand: "-" });
    assert.ok(hit);
    assert.equal(hit.strand, "-");
    assert.equal(hit.kind, "i-motif");
    assert.ok(hit.strandSequence.startsWith("CCC"));
  });

  it("returns nothing rather than the wrong molecule", () => {
    // A silent `members[0]` fallback is how the card and the picture came to
    // disagree. No match is a correct outcome; a plausible wrong one is not.
    const r = geometryOnly(G, "DNA", conditions());
    assert.equal(memberForCard(r.members, { kind: "i-motif", strand: "+" }), undefined);
    assert.equal(memberForCard(r.members, { kind: "triplex", strand: "+" }), undefined);
  });

  it("honours the topology within a strand", () => {
    const r = geometryOnly(G, "DNA", conditions());
    const anti = memberForCard(r.members, {
      kind: "g-quadruplex", strand: "+", topology: "antiparallel basket",
    });
    assert.equal(anti?.topology, "antiparallel");
  });

  it("every member carries the strand it was built from", () => {
    for (const seq of [G, C]) {
      for (const m of geometryOnly(seq, "DNA", conditions()).members) {
        assert.ok(m.strand === "+" || m.strand === "-", `${m.id} has no strand`);
        assert.ok(m.strandSequence.length > 0, `${m.id} has no strand sequence`);
        assert.ok(m.id.includes(`@${m.strand}`), `${m.id} id does not encode its strand`);
      }
    }
  });
});

describe("claim basis", () => {
  it("badges a genomic proxy differently from a measurement", () => {
    const proxy = { claim_basis: { target_semantics: "genomic_proxy", biophysically_grounded: false } };
    const measured = { claim_basis: { target_semantics: "biophysical", biophysically_grounded: true } };
    assert.equal(claimBadge(proxy as EvidenceCard).label, "genomic proxy");
    // Not the bare word "measured". The number on the card is a prediction
    // whose training *labels* were measurements, which is a claim about the
    // corpus and not about this value — and the atlas neighbours, which really
    // are measurements, need that word for themselves.
    assert.equal(
      claimBadge(measured as EvidenceCard).label,
      "prediction · measured training data",
    );
    assert.ok(!/^measured$/.test(claimBadge(measured as EvidenceCard).label));
    assert.notEqual(claimBadge(proxy as EvidenceCard).tone, claimBadge(measured as EvidenceCard).tone);
  });
});

describe("condition axis", () => {
  it("exposes all eight fields the atlas records", () => {
    assert.deepEqual(Object.keys(DEFAULT_CONDITIONS).sort(), [
      "crowder_pct", "k", "li_nh4", "mg", "na", "ph", "strand_conc", "temperature",
    ]);
  });

  it("keeps K+ and Na+ as separate fields that total correctly", () => {
    const c = conditions({ k: 100, na: 40, li_nh4: 10 });
    assert.equal(monovalent(c), 150);
    assert.equal(c.k, 100);
    assert.equal(c.na, 40);
  });
});

describe("out-of-domain values are withheld by default", () => {
  const card = (id: string, kind: string, inDomain: boolean): EvidenceCard =>
    ({
      id, kind, topology: "t", strand: "+", strand_sequence: "ACGT",
      source: "head:test",
      claim_basis: { target_semantics: "biophysical", biophysically_grounded: true },
      applicability: { in_domain: inDomain, warnings: inDomain ? [] : ["outside training range"] },
    }) as EvidenceCard;

  const cards = [card("g-quadruplex", "g-quadruplex", true), card("i-motif", "i-motif", false)];

  it("hides an extrapolation unless the user asks for it", () => {
    assert.deepEqual(visibleFolds(cards, false).map((c) => c.id), ["g-quadruplex"]);
    assert.deepEqual(extrapolatedFolds(cards).map((c) => c.id), ["i-motif"]);
  });

  it("shows everything once asked", () => {
    assert.equal(visibleFolds(cards, true).length, cards.length);
  });
});

describe("structure provenance", () => {
  it("names both solved topologies for the telomeric 22-mer", () => {
    const hits = templatesFor(HTELO);
    assert.equal(hits.length, 2, "one sequence, two solved folds chosen by cation");
    const topologies = hits.map((t) => t.topology);
    assert.ok(topologies.some((t) => /parallel propeller/.test(t)));
    assert.ok(topologies.some((t) => /antiparallel basket/.test(t)));
  });

  it("matches on the exact sequence only — a near match is a different molecule", () => {
    assert.equal(templatesFor(HTELO + "A").length, 0);
    assert.equal(templatesFor(HTELO.slice(1)).length, 0);
  });

  it("calls the drawing a schematic archetype whether or not a template exists", () => {
    assert.match(provenanceLabel(HTELO, "g-quadruplex"), /schematic structural archetype/);
    assert.match(provenanceLabel("ACGTACGTACGT", "g-quadruplex"), /schematic structural archetype/);
  });

  it("carries a verification note on every entry", () => {
    // Three of the five original entries were wrong — 2O3M was labelled c-MYC
    // when it is c-KIT, 2KYP was given c-KIT1's sequence when it is a 21-nt
    // c-kit2 structure, and 1EL2 was listed as an unmodified X-ray i-motif when
    // it is a 5-methylcytosine/uracil-substituted NMR structure. A wrong
    // accession beside a picture is a citation the reader will trust, so every
    // surviving entry states when it was checked against the deposition.
    assert.ok(TEMPLATES.length > 0);
    for (const t of TEMPLATES) {
      assert.ok(t.verified, `${t.pdb} has no verification note`);
      assert.match(t.verified, /rcsb\.org/, `${t.pdb} does not cite the deposition`);
    }
  });
});

describe("capability table is self-consistent", () => {
  it("gives every fold a reason, enabled or not", () => {
    for (const [kind, cap] of Object.entries(FOLD_CAPABILITIES)) {
      assert.ok(cap.reason.length > 30, `${kind} needs a stated reason`);
    }
  });

  it("still lists what it withheld, so absence is not read as 'not applicable'", () => {
    assert.equal(disabledFolds().length, 3);
  });
});

describe("the render path receives what it draws", () => {
  // The v0.4.3 regression, as a test. The store stopped producing a local
  // ensemble; ViewerPane still read `result.ensemble`, got `[]`, and the 3D
  // pane showed "No structure" — connected and offline alike — while the panels
  // beside it showed correct predictions. No type error fired, because an
  // optional array is a valid empty array.
  it("puts drawable structures where the viewer reads them", () => {
    const r = geometryOnly(HTELO, "DNA", conditions());
    assert.ok(r.members.length > 0, "members is what ViewerPane consumes");
    for (const m of r.members) {
      assert.ok(m.nucleotides.length > 0, `${m.id} has no coordinates to draw`);
      assert.ok(m.strands.length > 0, `${m.id} has no strand paths`);
    }
  });

  it("never ranks structures by a locally computed weight", () => {
    // NaScene filtered overlay members on `m.probability > 0.08`. That ranked
    // the pictures by a local Boltzmann weight — the exact thing the viewer must
    // not do — and once the field was removed the comparison became
    // `undefined > 0.08`, silently false, so overlay drew nothing.
    const r = geometryOnly(HTELO, "DNA", conditions());
    for (const m of r.members) {
      assert.ok(!("probability" in m), `${m.id} still carries a probability`);
    }
  });

  it("resolves a service card, using the topology it names", () => {
    const r = geometryOnly(HTELO, "DNA", conditions());
    const picked = memberForCard(r.members, {
      kind: "g-quadruplex", strand: "+", topology: "hybrid",
    });
    assert.ok(picked, "a service card must resolve to something drawable");
    assert.equal(picked.kind, "g-quadruplex");
    assert.equal(picked.topology, "hybrid");
  });

  it("refuses a G4 card that names no topology, because three archetypes fit", () => {
    // Ambiguity plus no match means no answer. `?? pool[0]` made "no match"
    // silently mean "parallel", since parallel is emitted first.
    const r = geometryOnly(HTELO, "DNA", conditions());
    assert.equal(memberForCard(r.members, { kind: "g-quadruplex", strand: "+" }), undefined);
    assert.equal(
      memberForCard(r.members, { kind: "g-quadruplex", strand: "+", topology: "chair" }),
      undefined,
    );
  });
});

describe("no ensemble language survives in the interface layer", () => {
  it("the local members describe themselves without stability claims", () => {
    const r = geometryOnly(HTELO, "DNA", conditions());
    for (const m of r.members) {
      assert.doesNotMatch(m.description, /probability|Boltzmann|reweight|ΔG|free energy/i);
      assert.match(m.description, /archetype/i);
    }
  });
});

describe("the preset library does not promise what the tool cannot do", () => {
  // It advertised "FMRP / DHX36 binders" for a sequence the service refuses,
  // and Mg2+ occupancy for an AC-motif no head models. A preset list is a set
  // of promises, and this checks each against what the tool actually produces.
  it("every preset's coverage badge matches its drawable structures", async () => {
    const { EXAMPLES } = await import("./examples.ts");
    for (const ex of EXAMPLES) {
      const r = geometryOnly(ex.sequence, ex.polymer, conditions());
      const modelled = r.members.some(
        (m) => m.kind === "g-quadruplex" || m.kind === "i-motif",
      );
      const coverage = ex.coverage ?? "full";

      if (coverage === "refused") {
        assert.ok(r.refusal, `${ex.name} is badged refused but was answered`);
        continue;
      }
      assert.ok(!r.refusal, `${ex.name} is not badged refused but was refused`);
      if (coverage === "none") {
        assert.equal(modelled, false,
          `${ex.name} is badged "not modelled" but draws a covered structure`);
      } else {
        assert.ok(modelled,
          `${ex.name} is badged "${coverage}" but draws no covered structure`);
      }
    }
  });

  it("no preset note claims protein binders or RNA analysis", async () => {
    const { EXAMPLES } = await import("./examples.ts");
    for (const ex of EXAMPLES) {
      assert.doesNotMatch(ex.note, /binder|FMRP|DHX36/i, ex.name);
    }
  });
});

describe("hybrid is a real archetype, not a relabelled parallel one", () => {
  // Human telomeric DNA in K+ adopts the hybrid (3+1) form, so this is the
  // common case, not an exotic one. It was missing entirely and a hybrid card
  // resolved to the parallel archetype — a different fold, drawn confidently.
  it("emits all three G4 topologies", () => {
    const r = geometryOnly(HTELO, "DNA", conditions());
    const topos = r.members.filter((m) => m.kind === "g-quadruplex").map((m) => m.topology);
    assert.deepEqual([...topos].sort(), ["antiparallel", "hybrid", "parallel"]);
  });

  it("builds hybrid with three columns one way and one the other", async () => {
    const { buildG4 } = await import("./geometry.ts");
    const runDirection = (topology: "parallel" | "antiparallel" | "hybrid") => {
      const g = buildG4(HTELO, topology);
      const tetrads = g.nucleotides.filter((n) => n.role === "tetrad");
      // Group the tetrad residues into their four columns, in sequence order,
      // and ask whether each column ascends or descends in y.
      const perColumn = 4;
      const size = Math.floor(tetrads.length / perColumn);
      const dirs: boolean[] = [];
      for (let c = 0; c < perColumn; c++) {
        const col = tetrads.slice(c * size, (c + 1) * size);
        if (col.length < 2) continue;
        dirs.push(col[col.length - 1]!.y > col[0]!.y);
      }
      return dirs;
    };

    const par = runDirection("parallel");
    const anti = runDirection("antiparallel");
    const hyb = runDirection("hybrid");

    assert.ok(par.every((up) => up), `parallel: all columns up, got ${par}`);
    assert.equal(anti.filter((up) => up).length, 2, `antiparallel: 2 up 2 down, got ${anti}`);
    assert.equal(hyb.filter((up) => up).length, 3, `hybrid is 3+1, got ${hyb}`);
    assert.notDeepEqual(hyb, par, "hybrid must differ from parallel");
    assert.notDeepEqual(hyb, anti, "hybrid must differ from antiparallel");
  });
});

describe("clicking an evidence card selects the molecule the card describes", () => {
  /**
   * The reviewer's exact reproduction. `EvidencePanels` stored the *service*
   * id (`i-motif`); local ids are `i-motif@-`; both viewers then fell back to
   * `members[0]`. Result:
   *
   *     click i-motif       -> g-quadruplex-parallel@+
   *     click g-quadruplex  -> g-quadruplex-parallel@+
   *
   * The strand category error the backend routing had just been fixed to
   * prevent, reintroduced one layer up by an id namespace mismatch. This is
   * that path, exercised through the same function the click handler calls.
   */
  const G = "AGGGTTAGGGTTAGGGTTAGGG";

  const serviceCards = [
    { kind: "g-quadruplex", strand: "+" as const, topology: "hybrid" },
    { kind: "i-motif", strand: "-" as const, topology: "intercalated C:C+ tetraplex" },
  ];

  it("each card resolves to its own molecule, not to whichever is first", () => {
    const r = geometryOnly(G, "DNA", conditions());
    const picked = serviceCards.map((c) => memberForCard(r.members, c));

    assert.equal(picked[0]?.id, "g-quadruplex-hybrid@+");
    assert.equal(picked[1]?.id, "i-motif@-");
    assert.notEqual(picked[0]?.id, picked[1]?.id, "two cards must not select one member");
  });

  it("the selected member carries the card's own strand and sequence", () => {
    const r = geometryOnly(G, "DNA", conditions());
    for (const card of serviceCards) {
      const m = memberForCard(r.members, card);
      assert.ok(m, `${card.kind} did not resolve`);
      assert.equal(m.strand, card.strand);
      assert.equal(m.kind, card.kind);
    }
    // and the i-motif really is the complement, not the typed strand
    const im = memberForCard(r.members, serviceCards[1]!)!;
    assert.ok(im.strandSequence.startsWith("CCC"), im.strandSequence);
    assert.equal(im.strandSequence.includes("GGG"), false);
  });

  it("an id that matches no member selects nothing at all", () => {
    // Both viewers now look up by exact id with no `?? members[0]`, so a stale
    // or foreign id shows an empty pane rather than a confident wrong picture.
    const r = geometryOnly(G, "DNA", conditions());
    assert.equal(r.members.find((m) => m.id === "i-motif"), undefined,
      "the bare service id must not accidentally match a local member");
    assert.equal(r.members.find((m) => m.id === ""), undefined);
  });
});

it("admits no per-card folding probability on the client type", async () => {
  // A branch of this project reintroduced `probabilityOfFolding` on every card,
  // added a hardcoded 0.40 for hairpin and 0.95 for duplex — both badged
  // `biophysically_grounded: true` — and divided the lot by their sum to make
  // "population percentages". That is the withdrawn ensemble in a cruder form:
  // not a Boltzmann weighting over a state space, just unlike quantities
  // normalised until they looked like one.
  //
  // The optional field on the client type is what let it back in without a
  // compile error, exactly as `ensemble: []` did a release earlier. This test
  // fails if the declaration returns.
  const src = await readFile(
    new URL("../quadcond/client.ts", import.meta.url),
    "utf8",
  );
  assert.ok(
    !/^\s*probabilityOfFolding\??:/m.test(src),
    "probabilityOfFolding is back on the card type",
  );
  assert.ok(
    !/^\s*population(Proportion|Percentage)\??:/m.test(src),
    "a normalised population share is back on the card type",
  );
});

describe("comparison exports", () => {
  it("writes a provenance header a spreadsheet and pandas both survive", async () => {
    // The header is commented so `pandas.read_csv(..., comment="#")` still
    // parses the table, and it is in the same file as the rows because a
    // separate provenance file is a file that gets left behind when the table
    // is emailed.
    const { batchToCsv } = await import("../quadcond/export.ts");
    const csv = batchToCsv(
      {
        quadcond_version: "0.4.7",
        scan: "batch",
        run_record: {
          workflow: "batch",
          quadcond_model_version: "0.4.8",
          model_artifact_sha256: "abc123",
          atlas_fingerprint: "fp123",
          prediction_schema_version: "1",
          input_mode: "fasta",
          n_sequences_in: 2,
          n_sequences_scored: 1,
          n_sequences_excluded: 1,
          n_rows: 1,
          conditions: [{ k: 140 }],
        },
        results: [
          {
            id: "tel22",
            record_index: 0,
            sequence: "AGGGTTAGGG",
            condition: "K140",
            predictions: {
              g4_tm: {
                value: 72.5,
                interval: [60, 85],
                target_semantics: "biophysical",
                applicability: { in_domain: true, warnings: [] },
              },
            },
          } as never,
        ],
        excluded: [
          { id: "junk", record_index: 1, reason: "no recognisable nucleotides" },
        ],
      },
      ["g4_tm"],
    );
    const header = csv.split("\n").filter((l) => l.startsWith("#"));
    assert.ok(header.some((l) => l.includes("model_version: 0.4.8")));
    // The artifact hash names the *file*; the version names a release. Every
    // export carries both now, from one shared run-record contract rather than
    // three hand-written header blocks that drifted apart.
    assert.ok(header.some((l) => l.includes("model_artifact_sha256: abc123")));
    assert.ok(header.some((l) => l.includes("atlas_fingerprint: fp123")));
    assert.ok(header.some((l) => l.includes("prediction_schema_version: 1")));
    // Exclusions travel in the file, not just in the panel — keyed by input
    // index as well as name, because identifiers are the user's labels and two
    // records may share one.
    assert.ok(
      header.some((l) => l.includes("junk") && l.includes("record_2")),
      "a batch that drops a sequence must say so in the export",
    );
    const body = csv.split("\n").filter((l) => l && !l.startsWith("#"));
    assert.equal(body[0], [
      "id", "record_index", "sequence", "condition", "head", "claim_basis",
      "value", "interval_low", "interval_high", "applicability", "notes",
    ].join(","));
    assert.ok(
      body[1]?.startsWith("tel22,0,AGGGTTAGGG,K140,g4_tm,biophysical,72.5,60,85,in_domain"),
    );
  });

  it("leaves an absent delta blank rather than zero", async () => {
    // A substitution that destroys the motif has no delta by construction.
    // Writing 0 would put "no effect" into a column a reader will average.
    const { mutationScanToCsv } = await import("../quadcond/export.ts");
    const csv = mutationScanToCsv({
      quadcond_version: "0.4.7",
      model_version: "0.4.7",
      scan: "mutation",
      wild_type: {
        sequence: "AGGGTTAGGG",
        complement: "CCCTAACCCT",
        length: 10,
        condition: "K140",
        condition_detail: { k: 140 },
        condition_imputed_fields: [],
        predictions: {},
        strand_by_head: { g4_tm: "+" },
        strand_sequence: { "+": "AGGGTTAGGG", "-": "CCCTAACCCT" },
      },
      heads: { g4_tm: { units: "degC" } as never },
      substitutions: [
        {
          position: 1, position_1based: 2, wild_type_base: "G", mutant_base: "A",
          label: "G2A", sequence: "AAGGTTAGGG", complement: "CCCTAACCTT",
          heads: {
            g4_tm: {
              strand: "+",
              motif: { state: "motif_lost", before: 1, after: 0 },
              wild_type: { state: "in_domain" },
              mutant: { state: "in_domain" },
            },
          },
        },
      ],
      ranking_note: "no significance test",
      delta_note: "not a validated interval on the difference",
      run_record: {
        workflow: "mutation_scan",
        quadcond_model_version: "0.4.8",
        model_artifact_sha256: "abc123",
        atlas_fingerprint: "fp123",
        prediction_schema_version: "1",
      },
    } as never);
    const lines = csv.split("\n");
    const header = lines.filter((l) => l.startsWith("#"));
    // The mutation export used to carry neither the artifact hash nor the
    // statement of what the paired spread is not — so a reader got a
    // `delta_paired_sd` column and no way to know it is not a validated
    // interval on the difference.
    assert.ok(header.some((l) => l.includes("model_artifact_sha256: abc123")));
    assert.ok(header.some((l) => l.includes("atlas_fingerprint: fp123")));
    assert.ok(
      header.some((l) => l.includes("not a validated interval")),
      "the export must carry the limitation, not just the number",
    );
    const row = lines.filter((l) => l && !l.startsWith("#"))[1]!;
    const cells = row.split(",");
    assert.equal(cells[7], "motif_lost");
    // wild_type_value, mutant_value, delta, delta_paired_sd, n_estimators
    assert.deepEqual(cells.slice(10), ["", "", "", "", ""]);
  });

  it("carries the strand and the provenance in a condition export", async () => {
    // A condition CSV whose rows do not say which strand they describe cannot
    // be read back without the panel that produced it — and the two series in
    // one sweep are, at a duplex position, two different molecules.
    const { conditionScanToCsv } = await import("../quadcond/export.ts");
    const csv = conditionScanToCsv({
      quadcond_version: "0.4.8",
      model_version: "0.4.8",
      scan: "condition",
      sequence: "GGGTTAGGGTTAGGGTTAGGG",
      complement: "CCCTAACCCTAACCCTAACCC",
      strand_sequence: { "+": "GGGTTAGGGTTAGGGTTAGGG", "-": "CCCTAACCCTAACCCTAACCC" },
      axis: "k",
      axis_label: "K+ (mM)",
      values: [50, 100],
      base_condition: { k: 100 },
      condition_imputed_fields: [],
      run_record: {
        workflow: "condition_scan",
        quadcond_model_version: "0.4.8",
        model_artifact_sha256: "abc123",
        atlas_fingerprint: "fp123",
        prediction_schema_version: "1",
      },
      summary_note: "in-domain points only",
      neighbour_tiers: ["experimental"],
      series: {
        g4_tm: {
          head: "g4_tm", kind: "G4", task: "regression", target: "tm",
          units: "degC", claim: "", target_semantics: "biophysical",
          biophysically_grounded: true,
          strand: "+", strand_sequence: "GGGTTAGGGTTAGGGTTAGGG",
          observed_response: {
            verdict: "responds", quantity: "value", spread: 12,
            error_scale: 3, summary_basis: "in_domain_points_only",
            n_points: 2, n_supported: 2, n_out_of_domain: 0,
          },
          points: [
            { value_of_axis: 50, strand: "+", value: 60,
              applicability: { state: "in_domain" } },
            { value_of_axis: 100, strand: "+", value: 72,
              applicability: { state: "in_domain" } },
          ],
        },
      },
      inert_heads: [],
      measured_neighbours: [],
      neighbour_note: "",
    } as never);
    const header = csv.split("\n").filter((l) => l.startsWith("#"));
    assert.ok(header.some((l) => l.includes("model_artifact_sha256: abc123")));
    assert.ok(header.some((l) => l.includes("atlas_fingerprint: fp123")));
    assert.ok(header.some((l) => l.includes("strand_by_head")));
    assert.ok(header.some((l) => l.includes("neighbour_tiers: experimental")));
    const body = csv.split("\n").filter((l) => l && !l.startsWith("#"));
    assert.equal(body[0], [
      "head", "units", "claim_basis", "strand", "strand_sequence", "k",
      "value", "interval_low", "interval_high", "folded_fraction",
      "applicability", "response_verdict", "response_spread", "head_error_scale",
    ].join(","));
    assert.ok(body[1]?.startsWith("g4_tm,degC,biophysical,+,GGGTTAGGGTTAGGGTTAGGG,50,60"));
  });

  it("names an export so it survives a rename", async () => {
    const { exportName } = await import("../quadcond/export.ts");
    const name = exportName("mutation-scan", "AGGGTTAGGGTTAGGG", "csv");
    assert.match(name, /^aenna_mutation-scan_AGGGTTAGGGTT_\d{4}-\d{2}-\d{2}\.csv$/);
  });
});

describe("the comparison panels", () => {
  const read = async (rel: string) =>
    readFile(new URL(rel, import.meta.url), "utf8");

  it("puts two axes on a variant without ever combining them", async () => {
    // A regulatory-effect prediction and a melting-temperature delta are model
    // output on different scales, and no weighting between them has been
    // calibrated. The panel shows both and ranks by the smaller of the two
    // percentile ranks; it must not grow a single blended number, which is the
    // shape of claim this project has already had to withdraw twice.
    const src = await read("../../components/panels/VariantPanel.tsx");
    for (const banned of ["combinedScore", "combined_score", "priorityScore",
                          "joint_score"]) {
      assert.ok(!src.includes(banned), `${banned} is back in the variant panel`);
    }
    assert.ok(/combined_rank/.test(src), "the rank the service computed is what is shown");

    const client = await read("../quadcond/client.ts");
    assert.ok(
      /regulatory:\s*\{[\s\S]{0,400}\} \| null/.test(client),
      "an absent regulatory record must be typed as null, not as a zero score",
    );
  });

  it("shows which model file answered and whether the release describes it", async () => {
    // A run record that reported `model_artifact_sha256: ""` named a version
    // and nothing that identifies the bytes. While the release artifact's
    // identity is disputed, the interface rendering its numbers must not be
    // silent about which artifact produced them.
    const src = await read("../../components/panels/VariantPanel.tsx");
    assert.ok(/variant-model-identity/.test(src), "model identity must be rendered");
    assert.ok(/model_file_sha256/.test(src), "identify the file by its contents");
    assert.ok(
      /matches_manifest === true/.test(src) && /matches_manifest === false/.test(src),
      "the three-valued manifest check must be rendered as three states",
    );
    assert.ok(
      /not a pass/.test(src),
      "an unavailable manifest comparison must not render as a pass",
    );
  });

  it("does not present a condition-inert head as condition-responsive", async () => {
    // im_pht is `fixed` on every axis: 160 training rows in a single buffer.
    // Rendering one overall "structural" verdict under a stated buffer implies
    // the buffer reached both heads.
    const src = await read("../../components/panels/VariantPanel.tsx");
    assert.ok(/variant-condition-responsiveness/.test(src));
    assert.ok(
      /no learned response to any condition axis/.test(src),
      "a wholly inert head must say so in those terms",
    );
    assert.ok(/reference-condition prediction/.test(src));
  });

  it("keeps motif-gain and motif-loss candidates visible without a delta", async () => {
    // These have no number by construction, so a delta-ranked table sorts them
    // to the bottom and off the end of a shortlist. That is a retrieval
    // failure, and the fix is a separate categorical list, not an invented
    // number.
    const src = await read("../../components/panels/VariantPanel.tsx");
    assert.ok(/variant-candidates/.test(src), "candidates must have their own view");
    assert.ok(/structural_candidates/.test(src));
    assert.ok(
      !/candidate[\s\S]{0,200}delta:/.test(src),
      "a categorical candidate must never be given a delta",
    );
  });

  it("says where the regulatory axis came from, or that there isn't one", async () => {
    // An empty regulatory column has two causes — the source found nothing, or
    // the source failed — and only one of them is about the locus.
    const src = await read("../../components/panels/VariantPanel.tsx");
    assert.ok(/variant-source/.test(src), "the source must be rendered");
    assert.ok(/regulatory_error/.test(src), "a failed source must be distinguishable");
    assert.ok(
      !/apiKey|api_key|ALPHAGENOME_API_KEY/.test(src),
      "the client must never hold or post an API key; the operator configures the source",
    );
  });

  it("invalidates a variant join when its window or coordinate changes", async () => {
    const src = await read("../../components/panels/VariantPanel.tsx");
    assert.ok(/const stale =/.test(src), "no staleness guard");
    assert.ok(/runId/.test(src), "no guard against an obsolete in-flight response");
    assert.ok(/locusKey/.test(src), "the coordinate must be part of what is pinned");
  });

  it("never calls the quarantined heuristic", async () => {
    // The branch panel this replaces imported `predict` from `./predict` and
    // ranked substitutions by its ΔΔG. That function is the one that cannot
    // tell K+ from Na+, whose calibration was never applied, and whose
    // reference state is not zero. `geometry-only.ts` exists so nothing in the
    // running application reaches it.
    for (const rel of [
      "../../components/panels/MutationExplorer.tsx",
      "../../components/panels/ConditionResponse.tsx",
      "../../components/panels/BatchPanel.tsx",
      "../../components/panels/Workbench.tsx",
    ]) {
      const src = await read(rel);
      assert.ok(
        !/from ["'].*na\/predict["']/.test(src),
        `${rel} imports the quarantined energy function`,
      );
      assert.ok(
        !/\bdeltaG\b|\bddG\b/.test(src),
        `${rel} mentions a free-energy difference; deltas come from the head that owns the quantity`,
      );
    }
  });

  it("breaks the condition line at unsupported points, not only at nulls", async () => {
    // The segmentation checked `p.y === null` alone, so an out-of-domain point
    // that still carried a number was joined into the path and merely given an
    // amber dot. A drawn line asserts a supported response between its
    // endpoints; a dot colour does not retract that.
    const src = await read("../../components/panels/ConditionResponse.tsx");
    assert.ok(
      /drawable\s*=\s*\(p[^)]*\)\s*=>\s*[\s\S]{0,120}p\.ok/.test(src),
      "line segmentation must consult applicability, not just the value",
    );
    assert.ok(
      /condition-show-unsupported/.test(src),
      "extrapolation must be an explicit, separate option",
    );
  });

  it("invalidates a condition sweep when its inputs change", async () => {
    // MutationExplorer has had this guard since v0.4.5. ConditionResponse had
    // none, so a completed sweep stayed on screen under a new sequence or
    // buffer, and a slow response could land after the user moved on and
    // repopulate the panel with the old scan.
    const src = await read("../../components/panels/ConditionResponse.tsx");
    assert.ok(/const stale =/.test(src), "no staleness guard");
    assert.ok(/runId/.test(src), "no guard against an obsolete in-flight response");
    assert.ok(
      /condition-frozen-context/.test(src),
      "the sweep must display the sequence and buffer it was run under",
    );
  });

  it("reads the scan-level delta note rather than a retired per-cell field", async () => {
    const src = await read("../../components/panels/MutationExplorer.tsx");
    assert.ok(
      /scan\.delta_note/.test(src),
      "the panel must read the note the service actually sends",
    );
    const stripped = src.replace(/\/\*[\s\S]*?\*\//g, "").replace(/^\s*\/\/.*$/gm, "");
    assert.ok(
      !/heads\[head\]\?\.delta_note/.test(stripped),
      "the per-cell delta_note was retired in v0.4.6 and always read undefined",
    );
  });

  it("calls nothing significant", async () => {
    // A fixed magnitude cutoff is a threshold with the word attached. There is
    // no null distribution here to test against, and the panel says so instead
    // of implying one.
    //
    // Comments are stripped first: this is a claim about what the interface
    // renders, not about what the source is allowed to discuss. The first
    // version failed on the docstring that explains why the label is gone,
    // which would have made the honest thing to do — writing down what was
    // wrong — the thing that breaks the build.
    const src = (await read("../../components/panels/MutationExplorer.tsx"))
      .replace(/\/\*[\s\S]*?\*\//g, "")
      .replace(/^\s*\/\/.*$/gm, "");
    assert.ok(
      !/significan/i.test(src),
      "no row may be labelled significant in the rendered panel",
    );
  });

  it("does not hardcode which heads exist", async () => {
    // A list of head names typed into a component goes stale the first time a
    // head is added and has no way to know it has — which is how the interface
    // came to show eleven heads while the model carried twelve.
    const src = await read("../../components/panels/Workbench.tsx");
    for (const name of ["g4_tm", "im_pht", "locus_peak_overlap_state", "g4_fold"]) {
      assert.ok(!src.includes(`"${name}"`), `Workbench hardcodes ${name}`);
    }
    assert.match(src, /meta\.workspaces/, "groups come from /info");
  });

  it("keeps the 3D canvas mounted when the Compare tab is shown", async () => {
    // Unmounting tore down the WebGL context on every switch: react-three-fiber
    // threw against a node Radix had already removed, and rebuilding the
    // context costs about four seconds on software rendering.
    const src = await read("../../components/layout/AppShell.tsx");
    assert.match(src, /forceMount[\s\S]{0,200}value="structure"/);
    assert.match(src, /data-\[state=inactive\]:hidden/);
  });
});
