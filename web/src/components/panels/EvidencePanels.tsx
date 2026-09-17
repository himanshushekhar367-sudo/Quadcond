/**
 * The panels that carry provenance: what the ensemble says, what was withheld,
 * and what the genomic proxies say in a place where they cannot be mistaken for
 * folding probabilities.
 *
 * The organising rule is that different kinds of statement get different real
 * estate. A melting temperature from a head trained on melting curves, an
 * antibody occupancy score from CUT&Tag, and a number from a local energy
 * function are three different claims, and an interface that stacks them in one
 * list has already merged them for the reader whatever the badges say.
 */
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import {
  allWarnings,
  claimBadge,
  extrapolatedFolds,
  visibleFolds,
  type EvidenceCard,
  type EvidenceResponse,
  type GenomicEvidenceItem,
  type SyntheticBanner,
} from "@/lib/quadcond/client";
import { disabledFolds } from "@/lib/na/capabilities";
import { memberForCard } from "@/lib/na/geometry-only";
import { MODEL_DISCLAIMER } from "@/lib/na/geometry-params";
import { pdbUrl, templatesForKind } from "@/lib/na/templates";
import type { StructureKind } from "@/lib/na/types";
import { useNA } from "@/lib/na/store";
import { AlertTriangle, Ban, ChevronRight, Dna, Eye, EyeOff, Microscope } from "lucide-react";

const TONE: Record<string, string> = {
  solid: "bg-emerald-500/15 text-emerald-700 border-emerald-500/30",
  warn: "bg-amber-500/15 text-amber-700 border-amber-500/30",
  muted: "bg-surface-3 text-fg-subtle border-border",
};

export function ClaimBadge({ fold }: { fold: EvidenceCard }) {
  const b = claimBadge(fold);
  return (
    <Badge variant="outline" className={`text-[9px] tracking-wide ${TONE[b.tone]}`}>
      {b.label}
    </Badge>
  );
}

/**
 * Strand-resolved evidence cards. **Not an ensemble, and deliberately not bars.**
 *
 * The previous panel drew each fold's Boltzmann weight as a percentage bar, and
 * bars in a stack are read as shares of one thing. They were not: the partition
 * function behind them contained a G-quadruplex, an i-motif on the *opposite*
 * strand, and a coil — while omitting the duplex, which is the state that
 * actually competes with both at physiological conditions, along with any strand
 * concentration or stoichiometry term. Percentages over an incomplete state
 * space are not uncertain, they are about something else.
 *
 * So each structure gets its own card with its own measured-scale quantity — a
 * melting temperature, a transitional pH — its interval, its strand and its
 * applicability. Nothing is normalised, and nothing sums to one.
 */
/**
 * The synthetic-provenance warning, rendered beside the numbers it is about.
 *
 * The workbench has shown this since v0.4.6, but `/evidence` was the one
 * numerical endpoint that did not carry the banner at all — so a predictor
 * flagged as synthetic produced a flagged `/predict` response and an unflagged
 * `/evidence` response containing the same numbers in cards, which is the view
 * the app opens on. A safeguard that covers five endpoints out of six covers
 * the wrong five.
 */
export function SyntheticNotice({ banner }: { banner?: SyntheticBanner }) {
  if (!banner) return null;
  return (
    <div
      className="mb-3 rounded-[var(--radius-md)] border border-red-500/50 bg-red-500/10 p-3"
      data-testid="synthetic-banner"
    >
      <p className="text-xs font-semibold text-red-700">{banner.warning}</p>
      <p className="mt-1 text-[0.7rem] text-red-700/80">{banner.note}</p>
    </div>
  );
}

export function EvidenceCards({
  cards,
  note,
  syntheticModel,
}: {
  cards: EvidenceCard[];
  note?: string;
  syntheticModel?: SyntheticBanner;
}) {
  const selectedId = useNA((s) => s.selectedId);
  const selectStructure = useNA((s) => s.selectStructure);
  const show = useNA((s) => s.showExtrapolation);
  const setShow = useNA((s) => s.setShowExtrapolation);
  const members = useNA((s) => s.result.members);

  const shown = visibleFolds(cards, show);
  const hidden = extrapolatedFolds(cards);

  /**
   * A card click selects the *drawable member* for that card, never the card's
   * own id.
   *
   * Service ids are `i-motif`; local ids are `i-motif@-`. Storing the service id
   * meant nothing matched and both viewers fell back to `members[0]`, so
   * clicking the i-motif card drew a parallel G4 — the same strand category
   * error the backend routing had just been fixed to prevent, reintroduced one
   * layer up by an id namespace mismatch.
   */
  const pick = (card: EvidenceCard) => {
    const member = memberForCard(members, card);
    if (member) selectStructure(member.id);
  };
  const drawableFor = (card: EvidenceCard) => memberForCard(members, card);

  return (
    <div className="space-y-3">
      {/* Above the numbers, not below them. */}
      <SyntheticNotice banner={syntheticModel} />
      {note ? (
        <div className="rounded-lg border border-border bg-surface p-3 text-xs leading-relaxed text-fg-muted">
          <p>Separate predictions for each strand, not a thermodynamic ensemble. Values do not sum to one; duplex competition is not modelled.</p>
          <details className="mt-1">
            <summary className="text-fg-subtle">Interpretation details</summary>
            <p className="mt-2">{note}</p>
          </details>
        </div>
      ) : null}

      <ul className="space-y-3">
        {shown.map((c) => {
          const drawable = drawableFor(c);
          const active = drawable != null && drawable.id === selectedId;
          const out = c.applicability?.in_domain === false;
          const value =
            c.predictedTm != null
              ? { label: "melting temperature", text: `${c.predictedTm} °C`, ci: c.predictedTm_interval }
              : c.predictedPhT != null
                ? { label: "transitional pH", text: `pH ${c.predictedPhT}`, ci: c.predictedPhT_interval }
                : null;
          return (
            <li key={c.id}>
              <button
                type="button"
                data-testid="evidence-card"
                data-card-id={c.id}
                data-kind={c.kind}
                data-service-topology={c.topology}
                data-strand={c.strand}
                data-topology={drawable?.topology ?? ""}
                data-member-id={drawable?.id ?? ""}
                aria-pressed={active}
                onClick={() => pick(c)}
                disabled={!drawable}
                title={drawable ? undefined : "no archetype available for this topology"}
                className={`w-full rounded-xl border p-4 text-left transition-all duration-300 ${
                  active
                    ? "border-accent/50 bg-accent/5"
                    : "border-border bg-surface-2 hover:bg-surface-2 hover:border-border"
                } ${out ? "border-dashed opacity-80" : ""}`}
              >
                <div className="mb-3 flex items-center justify-between gap-2">
                  <div className="flex items-center gap-2">
                    <span className={`text-sm font-bold ${active ? "text-accent" : "text-fg"}`}>
                      {c.kind === "g-quadruplex" ? "G-quadruplex" : "i-motif"}
                    </span>
                    <ClaimBadge fold={c} />
                  </div>
                  {/* Which strand, always. A G4 and an i-motif at one duplex
                      position sit on opposite strands, and the previous build
                      evaluated both heads on whichever strand the user typed. */}
                  <Badge variant="outline" className="border-border font-mono text-[10px] text-fg-subtle">
                    strand {c.strand}
                  </Badge>
                </div>

                {value ? (
                  <div className="mb-2">
                    <p className="text-[10px] uppercase tracking-widest text-fg-subtle">{value.label}</p>
                    <p className="font-mono text-lg font-bold tabular-nums text-accent">
                      {value.text}
                      {value.ci ? (
                        <span className="ml-2 font-normal text-[11px] text-fg-subtle">
                          [{value.ci[0]}, {value.ci[1]}]
                        </span>
                      ) : null}
                    </p>
                  </div>
                ) : null}

                {c.foldedFraction != null ? (
                  <p className="text-[11px] leading-relaxed text-fg-muted">
                    Folded fraction of the <em>isolated strand</em>:{" "}
                    <span className="font-mono tabular-nums">{c.foldedFraction.toFixed(3)}</span>
                    <span className="block text-[10px] text-fg-subtle">{c.foldedFraction_note}</span>
                  </p>
                ) : null}

                {/* The prediction stands whether or not this viewer can draw
                    the topology it names. Saying so beats drawing a different
                    fold, and beats a card that does nothing when clicked. */}
                {!drawable ? (
                  <p className="mt-2 rounded border border-amber-500/25 bg-amber-500/5 px-2 py-1 text-[10px] leading-relaxed text-amber-700/90">
                    No 3D archetype for the <strong>{c.topology}</strong> topology, so
                    nothing is drawn for this card. The prediction above is unaffected.
                  </p>
                ) : null}

                <p className="mt-2 break-all border-l-2 border-border pl-2 font-mono text-[9px] leading-snug text-fg-subtle">
                  {c.strand_sequence}
                </p>
                <p className="mt-1 pl-2 font-mono text-[9px] text-fg-subtle">{c.source}</p>
              </button>
            </li>
          );
        })}
      </ul>

      {hidden.length > 0 ? (
        <div className="rounded-lg border border-amber-500/25 bg-amber-500/5 p-3">
          <div className="flex items-start gap-2">
            <AlertTriangle className="mt-0.5 h-3.5 w-3.5 shrink-0 text-amber-700" />
            <div className="space-y-2">
              <p className="text-[11px] leading-relaxed text-amber-700/90">
                {hidden.length} structure{hidden.length > 1 ? "s are" : " is"} hidden
                because the query falls outside what the head was trained on.
              </p>
              <ul className="space-y-1">
                {allWarnings(hidden).map((w) => (
                  <li key={w} className="text-[10px] leading-relaxed text-amber-700/70">
                    · {w}
                  </li>
                ))}
              </ul>
              <Button
                size="sm"
                variant="ghost"
                className="h-6 px-2 text-[10px] tracking-wider text-amber-700 hover:bg-amber-500/10"
                onClick={() => setShow(!show)}
              >
                {show ? <EyeOff className="mr-1 h-3 w-3" /> : <Eye className="mr-1 h-3 w-3" />}
                {show ? "Hide extrapolation" : "Show extrapolation anyway"}
              </Button>
            </div>
          </div>
        </div>
      ) : null}

      {shown.length > 0 && allWarnings(shown).length > 0 ? (
        <ul className="space-y-1 rounded-lg border border-border bg-surface-2 p-3">
          {allWarnings(shown).map((w) => (
            <li key={w} className="text-[10px] leading-relaxed text-fg-subtle">
              · {w}
            </li>
          ))}
        </ul>
      ) : null}
    </div>
  );
}

/**
 * Antibody occupancy, in its own panel.
 *
 * These scores are an experimental observation of a *different quantity* than
 * folding: whether BG4 or iMab bound a locus in live HEK293T cells. No buffer
 * was measured, the iMab antibody's specificity is contested in the literature,
 * and chromatin accessibility contributes to targeted CUT&Tag signal. None of
 * that stops the number from looking exactly like a probability of folding if
 * it is printed beside one, which is why it is not.
 */
export function GenomicEvidencePanel({
  items,
  note,
}: {
  items: GenomicEvidenceItem[];
  note?: string;
}) {
  const show = useNA((s) => s.showExtrapolation);
  const setShow = useNA((s) => s.setShowExtrapolation);
  if (!items.length) return null;

  const outOfDomain = items.filter((g) => !g.refused && g.applicability?.in_domain === false);
  // Refusals remain visible as explanations, with no score, under either toggle state.
  const visible = items.filter((g) => g.refused || show || g.applicability?.in_domain !== false);

  return (
    <details className="group [&_summary::-webkit-details-marker]:hidden">
      <summary className="flex cursor-pointer select-none items-center justify-between text-fg-subtle transition-colors hover:text-fg">
        <div className="flex items-center gap-2">
          <Microscope className="h-4 w-4 text-amber-700" />
          <p className="text-xs font-bold uppercase tracking-widest">Genomic evidence</p>
        </div>
        <ChevronRight className="h-4 w-4 transition-transform group-open:rotate-90" />
      </summary>

      <p className="mt-2 text-[11px] leading-relaxed text-amber-700/80">
        {note ??
          "Antibody occupancy at genomic loci. An observation of a different quantity than folding — never read as P(folds)."}
      </p>

      <ul className="mt-3 space-y-2">
        {visible.map((g) => (
          <li key={g.head} data-testid="genomic-readout" data-refused={g.refused === true} className={`rounded-lg border bg-surface-2 p-3 ${
            g.applicability?.in_domain === false
              ? "border-dashed border-amber-500/30"
              : "border-amber-500/20"
          }`}>
            <div className="flex items-center justify-between gap-2">
              <span className="text-xs font-semibold text-fg">{g.label}</span>
              <div className="flex items-center gap-2">
                <span className="font-mono text-[10px] text-fg-subtle">strand {g.strand}</span>
                <Badge variant="outline" className={`font-mono text-[10px] ${TONE.warn}`}>
                  {g.refused ? "Refused" : g.score.toFixed(3)}
                </Badge>
              </div>
            </div>
            <p className="mt-2 text-[10px] leading-relaxed text-fg-subtle">{g.claim}</p>
            {g.refused ? (
              <p className="mt-2 text-xs leading-relaxed text-amber-700" role="status">
                {g.refusal_reason} No score is available; extrapolation cannot override this refusal.
              </p>
            ) : g.applicability?.in_domain === false ? (
              <p className="mt-1 text-[10px] leading-relaxed text-amber-700/80">
                Extrapolation: {g.applicability.warnings[0]}
              </p>
            ) : null}
          </li>
        ))}
      </ul>

      {/*
        The same gate as everywhere else. These two scores were printed
        unconditionally, and at the app's own default buffer both heads are out
        of domain -- their training rows carry Mg2+ 0.5 / pH 7.2 against a UI
        default of Mg2+ 1 / pH 7.4 -- so the panel a user sees on first load was
        showing extrapolations with no toggle in front of them.
      */}
      {outOfDomain.length && !show ? (
        <div className="mt-3 rounded-lg border border-amber-500/25 bg-amber-500/5 p-3">
          <div className="flex items-start gap-2">
            <AlertTriangle className="mt-0.5 h-3.5 w-3.5 shrink-0 text-amber-700" />
            <div className="space-y-2">
              <p className="text-[11px] leading-relaxed text-amber-700/90">
                {outOfDomain.length} score{outOfDomain.length > 1 ? "s are" : " is"} hidden:
                the query falls outside what {outOfDomain.length > 1 ? "these heads were" : "this head was"} trained on.
              </p>
              <ul className="space-y-1">
                {[...new Set(outOfDomain.flatMap((g) => g.applicability?.warnings ?? []))].map((w) => (
                  <li key={w} className="text-[10px] leading-relaxed text-amber-700/70">
                    · {w}
                  </li>
                ))}
              </ul>
              <Button
                size="sm"
                variant="ghost"
                className="h-6 px-2 text-[10px] tracking-wider text-amber-700 hover:bg-amber-500/10"
                onClick={() => setShow(true)}
              >
                <Eye className="mr-1 h-3 w-3" />
                Show extrapolation anyway
              </Button>
            </div>
          </div>
        </div>
      ) : null}
    </details>
  );
}

/**
 * The folds this build refuses to draw, listed rather than silently absent.
 *
 * A missing panel reads as "not applicable to this sequence". A listed one with
 * a reason reads as "we do not know", which is the true statement.
 */
export function WithheldFolds() {
  const withheld = disabledFolds();
  if (!withheld.length) return null;
  return (
    <details className="group [&_summary::-webkit-details-marker]:hidden">
      <summary className="flex cursor-pointer select-none items-center justify-between text-fg-subtle transition-colors hover:text-fg">
        <div className="flex items-center gap-2">
          <Ban className="h-4 w-4 text-fg-subtle" />
          <p className="text-xs font-bold uppercase tracking-widest">
            Not shown ({withheld.length})
          </p>
        </div>
        <ChevronRight className="h-4 w-4 transition-transform group-open:rotate-90" />
      </summary>
      <ul className="mt-3 space-y-2">
        {withheld.map((w) => (
          <li key={w.kind} className="rounded-lg border border-border bg-surface-2 p-3">
            <p className="text-xs font-semibold text-fg-muted">{w.kind}</p>
            <p className="mt-1 text-[10px] leading-relaxed text-fg-subtle">{w.reason}</p>
          </li>
        ))}
      </ul>
    </details>
  );
}

/** RNA and anything else the tool declines outright. */
export function RefusalNotice({ reason }: { reason: string }) {
  return (
    <div className="rounded-xl border border-amber-500/40 bg-amber-500/10 p-4">
      <div className="flex items-start gap-2">
        <Dna className="mt-0.5 h-4 w-4 shrink-0 text-amber-700" />
        <div>
          <p className="text-xs font-bold uppercase tracking-widest text-amber-700">
            No prediction for this input
          </p>
          <p className="mt-2 text-[11px] leading-relaxed text-amber-700/80">{reason}</p>
        </div>
      </div>
    </div>
  );
}

/**
 * What the picture is, stated under the picture.
 *
 * "Schematic structural archetype" is the accurate name for what the geometry
 * builder produces: correct helical parameters, correct tetrad stacking,
 * correct intercalation spacing -- and nothing sequence-specific beyond where
 * the tracts fall. Loop conformation, groove widths and syn/anti glycosidic
 * angles are not predicted, and those are most of what distinguishes one real
 * quadruplex from another.
 *
 * Where a solved structure exists for the exact sequence, it is named and
 * linked rather than substituted. The user can then compare, which is more
 * useful than being shown a deposited structure and left unsure which parts of
 * the image came from an experiment.
 */
export function StructureProvenance({
  sequence,
  kind,
}: {
  sequence: string;
  kind: StructureKind;
}) {
  const matches = templatesForKind(sequence, kind);
  return (
    <div className="rounded-xl border border-border bg-surface-2 p-4">
      <p className="text-[10px] font-bold uppercase tracking-widest text-fg-subtle">
        What you are looking at
      </p>
      <p className="mt-2 text-[11px] leading-relaxed text-fg-muted">
        Schematic structural archetype. {MODEL_DISCLAIMER}
      </p>

      {matches.length ? (
        <div className="mt-3 border-t border-border pt-3">
          <p className="text-[10px] leading-relaxed text-emerald-700/90">
            {matches.length === 1
              ? "A solved structure exists for this exact sequence:"
              : `${matches.length} solved structures exist for this exact sequence — the fold depends on the buffer:`}
          </p>
          <ul className="mt-2 space-y-2">
            {matches.map((t) => (
              <li key={t.pdb} className="flex items-start justify-between gap-3">
                <div>
                  <a
                    href={pdbUrl(t.pdb)}
                    target="_blank"
                    rel="noreferrer noopener"
                    className="font-mono text-[11px] font-bold text-accent hover:underline"
                  >
                    {t.pdb}
                  </a>
                  <p className="text-[10px] leading-snug text-fg-subtle">{t.name}</p>
                </div>
                <span className="shrink-0 text-right text-[10px] text-fg-subtle">
                  {t.topology}
                  <br />
                  {t.method} · {t.conditions}
                </span>
              </li>
            ))}
          </ul>
        </div>
      ) : null}
    </div>
  );
}

/**
 * Peak-overlap posterior, gated like everything else.
 *
 * The panel this replaces was titled "Joint occupancy (locus_state)" and
 * rendered its four numbers unconditionally. Both halves were wrong. The head
 * does not measure joint occupancy — BG4 and iMab CUT&Tag were parallel
 * reactions on separate aliquots, so a "both" call is population-level peak
 * co-localisation and nothing in the design observes two structures on one
 * molecule. And every training window is 201 nt, so for the 22-mer a user
 * typically types, the posterior on screen was pure extrapolation with no
 * indication of it: the service was sending an applicability block and the
 * client type dropped it.
 */
export function LocusOverlapPanel({
  joint,
}: {
  joint: NonNullable<EvidenceResponse["locus_peak_overlap_state"]>;
}) {
  const show = useNA((s) => s.showExtrapolation);
  const setShow = useNA((s) => s.setShowExtrapolation);
  const out = joint.applicability?.in_domain === false;

  return (
    <details className="group [&_summary::-webkit-details-marker]:hidden">
      <summary className="flex cursor-pointer select-none items-center justify-between text-fg-subtle transition-colors hover:text-fg">
        <div className="flex items-center gap-2">
          <Microscope className="h-4 w-4 text-amber-700" />
          <p className="text-xs font-bold uppercase tracking-widest">
            Peak overlap (201-nt windows)
          </p>
        </div>
        <ChevronRight className="h-4 w-4 transition-transform group-open:rotate-90" />
      </summary>

      <div className="mt-3 rounded-xl border border-amber-500/20 bg-surface-2 p-4">
        <div className="mb-2 flex items-center justify-between gap-2">
          <Badge variant="outline" className={`text-[9px] ${TONE.warn}`}>
            genomic proxy
          </Badge>
          <span className="font-mono text-[10px] text-fg-subtle">
            strand {joint.strand}
          </span>
        </div>

        {joint.refused ? (
          <div data-testid="locus-refusal" role="status" className="space-y-2 text-xs leading-relaxed text-amber-700">
            <p className="font-semibold">Prediction refused</p>
            <p>{joint.refusal_reason}</p>
            <p>No posterior is available; extrapolation cannot override this refusal.</p>
          </div>
        ) : out && !show ? (
          <>
            <div className="flex items-start gap-2">
              <AlertTriangle className="mt-0.5 h-3.5 w-3.5 shrink-0 text-amber-700" />
              <div className="space-y-2">
                <p className="text-[11px] leading-relaxed text-amber-700/90">
                  Not shown: this query is outside what the head was trained on.
                </p>
                <ul className="space-y-1">
                  {(joint.applicability?.warnings ?? []).map((w) => (
                    <li key={w} className="text-[10px] leading-relaxed text-amber-700/70">
                      · {w}
                    </li>
                  ))}
                </ul>
                <Button
                  size="sm"
                  variant="ghost"
                  className="h-6 px-2 text-[10px] tracking-wider text-amber-700 hover:bg-amber-500/10"
                  onClick={() => setShow(true)}
                >
                  <Eye className="mr-1 h-3 w-3" />
                  Show extrapolation anyway
                </Button>
              </div>
            </div>
          </>
        ) : (
          <>
            {out ? (
              <p className="mb-2 text-[10px] leading-relaxed text-amber-700/80">
                Extrapolation — shown because you asked.{" "}
                {(joint.applicability?.warnings ?? [])[0]}
              </p>
            ) : null}
            <div className="space-y-1">
              {Object.entries(joint.posterior).map(([k, v]) => (
                <div key={k} className="flex items-center justify-between gap-3">
                  <span className="text-[11px] text-fg-muted">{k}</span>
                  <span className="font-mono text-[11px] tabular-nums text-accent">
                    {(v as number).toFixed(3)}
                  </span>
                </div>
              ))}
            </div>
          </>
        )}

        <p className="mt-3 border-t border-border pt-3 text-[10px] leading-relaxed text-fg-subtle">
          {joint.note}
        </p>
      </div>
    </details>
  );
}
