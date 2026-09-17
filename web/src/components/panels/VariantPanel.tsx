/**
 * Which substitution is worth measuring — asked of two models at once.
 *
 * A regulatory-effect prediction says a promoter variant should change
 * expression. It does not say how. A G-quadruplex melting-temperature delta
 * says the same edit should cost the local structure fifteen degrees. It does
 * not say whether anything downstream cares. The panel puts both on the same
 * substitution and lets the reader see which variants score on both axes.
 *
 * What it will not do is add them together. The two numbers are model output on
 * different scales, nothing has calibrated a weighting between them, and a
 * single combined score would be the third time this project has had to
 * withdraw a number that looked like a measurement. So the display is the two
 * axes, the quadrant against stated thresholds, and a rank that is the smaller
 * of the two — high only when both are high.
 *
 * The regulatory source is whatever the operator configured on the service.
 * The panel never asks for a file path or an API key, and when no source is
 * configured it says the axis is absent rather than drawing an empty column
 * that reads as "no effect".
 */
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { AlertTriangle, Crosshair, Download, Loader2 } from "lucide-react";
import { useNA } from "@/lib/na/store";
import {
  QuadcondUnavailable,
  scanVariant,
  type VariantRow,
  type VariantScan,
} from "@/lib/quadcond/client";
import { download, exportName, toJson, variantScanToCsv } from "@/lib/quadcond/export";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";

const QUADRANT_LABEL: Record<VariantRow["quadrant"], string> = {
  both: "regulatory + structural",
  structural_only: "structural only",
  regulatory_only: "regulatory only",
  neither: "neither",
  unclassified: "not classifiable",
};

const QUADRANT_TONE: Record<VariantRow["quadrant"], string> = {
  both: "border-emerald-500/40 bg-emerald-500/10 text-emerald-700",
  structural_only: "border-accent/40 bg-accent/10 text-fg",
  regulatory_only: "border-sky-500/40 bg-sky-500/10 text-sky-700",
  neither: "border-border bg-bg/40 text-fg-muted",
  unclassified: "border-dashed border-amber-500/40 bg-amber-500/5 text-amber-700",
};

/** The 2×2, drawn from the counts the service computed rather than recounted here. */
function Quadrants({ counts }: { counts: Record<string, number> }) {
  const cell = (key: VariantRow["quadrant"], title: string, sub: string) => (
    <div
      className={`rounded-[var(--radius-md)] border p-3 ${QUADRANT_TONE[key]}`}
      data-testid="variant-quadrant"
      data-quadrant={key}
    >
      <div className="text-lg font-semibold tabular-nums">{counts[key] ?? 0}</div>
      <div className="text-[0.7rem] font-medium">{title}</div>
      <div className="mt-0.5 text-[0.65rem] opacity-80">{sub}</div>
    </div>
  );
  return (
    <div className="mt-3 grid grid-cols-2 gap-2">
      {cell("both", "Both axes high", "candidate structure-mediated mechanism")}
      {cell("regulatory_only", "Regulatory only", "effect without a structural account here")}
      {cell("structural_only", "Structural only", "structure moves; no regulatory signal")}
      {cell("neither", "Neither", "low on both axes")}
      {cell("unclassified", "Not classifiable", "one or both numerical axes unavailable")}
    </div>
  );
}

export function VariantPanel() {
  const sequence = useNA((s) => s.sequence);
  const conditions = useNA((s) => s.conditions);
  const backend = useNA((s) => s.backend);
  const setSequence = useNA((s) => s.setSequence);
  const setPolymer = useNA((s) => s.setPolymer);

  const [chromosome, setChromosome] = useState("");
  const [start, setStart] = useState("");
  const [scan, setScan] = useState<VariantScan | null>(null);
  const [busy, setBusy] = useState(false);
  const [showAll, setShowAll] = useState(false);
  const [error, setError] = useState<string | null>(null);

  /** Same pin-and-discard guard as the other comparison panels. A join is only
   *  meaningful for the window it was computed on, and a coordinate the user
   *  has since edited is a different window. */
  const conditionKey = JSON.stringify(conditions);
  const [pin, setPin] = useState<{ sequence: string; conditions: string; locus: string } | null>(null);
  const locusKey = `${chromosome}:${start}`;
  const stale =
    pin !== null &&
    (pin.sequence !== sequence || pin.conditions !== conditionKey || pin.locus !== locusKey);

  useEffect(() => {
    if (stale) {
      setScan(null);
      setPin(null);
      setError(null);
    }
  }, [stale]);

  const runId = useRef(0);

  const run = useCallback(async () => {
    const id = ++runId.current;
    const pos = Number(start);
    if (!/^\d+$/.test(start) || !Number.isSafeInteger(pos) || pos < 1) {
      setError("Start must be the 1-based coordinate of the window's first base.");
      return;
    }
    const pinned = { sequence, conditions: conditionKey, locus: locusKey };
    setBusy(true);
    setError(null);
    try {
      const result = await scanVariant(sequence, { chromosome, start: pos }, conditions);
      if (id !== runId.current) return;
      setScan(result);
      setPin(pinned);
    } catch (err) {
      if (id !== runId.current) return;
      setError(err instanceof QuadcondUnavailable ? err.message : String(err));
      setScan(null);
    } finally {
      if (id === runId.current) setBusy(false);
    }
  }, [sequence, chromosome, start, conditions, conditionKey, locusKey]);

  const top = useMemo(
    () => (scan ? (showAll ? scan.variants : scan.variants.slice(0, 25)) : []),
    [scan, showAll],
  );

  if (backend !== "connected") return null;

  const source = scan?.run_record.regulatory_source;
  const identity = scan?.run_record.provenance?.model ?? null;
  const responsiveness = scan?.run_record.condition_responsiveness ?? null;
  const candidateTotal = scan
    ? Object.values(scan.structural_candidate_counts ?? {}).reduce((a, b) => a + b, 0)
    : 0;
  const units = scan && scan.run_record.structural_axis_head
    ? (scan.heads[scan.run_record.structural_axis_head]?.units ?? "")
    : "";

  return (
    <section
      className="rounded-[var(--radius-lg)] border border-border bg-surface p-4"
      data-testid="variant-panel"
    >
      <header className="flex flex-wrap items-center justify-between gap-2">
        <div>
          <h3 className="flex items-center gap-2 text-sm font-medium text-fg">
            <Crosshair className="h-4 w-4" aria-hidden />
            Variant join
          </h3>
          <p className="text-xs text-fg-subtle">
            Every substitution in this window, against its AlphaGenome regulatory
            score for the same edit.
          </p>
        </div>
        <div className="flex items-center gap-1">
          <Button size="sm" onClick={run} disabled={busy || !chromosome.trim() || !start.trim()} data-testid="run-variant-scan">
            {busy ? <Loader2 className="h-3.5 w-3.5 animate-spin" aria-hidden /> : null}
            {scan ? "Re-run" : "Join"}
          </Button>
          {scan ? (
            <>
              <Button
                size="sm"
                variant="ghost"
                onClick={() =>
                  download(exportName("variant-join", scan.window.sequence, "csv"),
                           variantScanToCsv(scan))
                }
              >
                <Download className="h-3.5 w-3.5" aria-hidden /> CSV
              </Button>
              <Button
                size="sm"
                variant="ghost"
                onClick={() =>
                  download(exportName("variant-join", scan.window.sequence, "json"),
                           toJson(scan), "application/json")
                }
              >
                JSON
              </Button>
            </>
          ) : null}
        </div>
      </header>
      <div className="mt-4 flex flex-wrap items-center gap-2 rounded-xl bg-surface-2 p-3">
        <Button size="sm" variant="outline" disabled={busy} onClick={() => {
          setPolymer("DNA");
          setSequence("AGGAGGGCAGAGAGCTGGGGCCTCGGACTCACCCGACGCTTGTGATGAGCTGCACCCAGGA");
          setChromosome("chr22");
          setStart("36201668");
          setScan(null); setPin(null); setShowAll(false); setError(null);
        }}>Use chr22 example</Button>
        <p className="text-xs text-fg-muted">Negative control · requires the matching saved Atlas export.</p>
      </div>

      <div className="mt-3 flex flex-wrap items-end gap-2">
        <label className="text-[0.7rem] text-fg-subtle">
          Chromosome
          <input
            placeholder="e.g. chr22"
            value={chromosome}
            onChange={(e) => setChromosome(e.target.value)}
            data-testid="variant-chromosome"
            className="mt-1 block w-28 rounded-[var(--radius-md)] border border-border bg-bg/60 px-2 py-1 font-mono text-xs text-fg focus:border-accent focus:outline-none"
          />
        </label>
        <label className="text-[0.7rem] text-fg-subtle">
          Start of window (1-based)
          <input
            placeholder="e.g. 36201668"
            value={start}
            onChange={(e) => setStart(e.target.value)}
            inputMode="numeric"
            data-testid="variant-start"
            className="mt-1 block w-40 rounded-[var(--radius-md)] border border-border bg-bg/60 px-2 py-1 font-mono text-xs text-fg focus:border-accent focus:outline-none"
          />
        </label>
        <p className="text-[0.65rem] text-fg-subtle">
          Use the verified GRCh38 forward-strand coordinate of the first base.
          This panel does not fetch the reference assembly to verify your input.
        </p>
      </div>

      {error ? (
        <p className="mt-3 rounded-[var(--radius-md)] border border-amber-500/40 bg-amber-500/10 p-2 text-xs text-amber-700">
          {error}
        </p>
      ) : null}

      {scan ? (
        <>
          {/* Where the regulatory axis came from, or that there isn't one. */}
          <p
            className="mt-3 rounded-[var(--radius-md)] border border-border bg-bg/40 p-2 text-[0.7rem] leading-relaxed text-fg-subtle"
            data-testid="variant-source"
          >
            <span className="font-mono text-fg">{source?.source ?? "unknown"}</span>{" "}
            {source?.note}
            {scan.run_record.regulatory_error ? (
              <span className="text-amber-700">
                {" "}
                The source failed for this request: {scan.run_record.regulatory_error}. An
                empty regulatory column below is that failure, not a quiet locus.
              </span>
            ) : null}
          </p>

          {/* Which model answered, by its contents, and whether the release
              describes it. An artifact whose identity is disputed must not be
              invisible in the interface that renders its numbers. */}
          {identity ? (
            <p
              className={`mt-2 rounded-[var(--radius-md)] border p-2 text-[0.7rem] leading-relaxed ${
                identity.matches_manifest === true
                  ? "border-border bg-bg/40 text-fg-subtle"
                  : "border-amber-500/40 bg-amber-500/10 text-amber-700"
              }`}
              data-testid="variant-model-identity"
            >
              Model{" "}
              <span className="font-mono text-fg">
                {identity.model_version_reported_by_artifact ?? "unknown"}
              </span>
              {identity.model_file_sha256 ? (
                <>
                  {" "}
                  · sha256{" "}
                  <span className="font-mono">
                    {identity.model_file_sha256.slice(0, 12)}…
                  </span>
                </>
              ) : null}{" "}
              ·{" "}
              {identity.matches_manifest === true
                ? "matches the release manifest"
                : identity.matches_manifest === false
                  ? "does NOT match the release manifest"
                  : "manifest comparison unavailable — not a pass"}
              . {identity.note}
            </p>
          ) : null}

          {/* Per head, per axis. One overall "structural" verdict hides that a
              head may have no learned response to the buffer that was set. */}
          {responsiveness ? (
            <details
              className="mt-2 rounded-[var(--radius-md)] border border-border bg-bg/40 p-2 text-[0.7rem]"
              data-testid="variant-condition-responsiveness"
            >
              <summary className="cursor-pointer text-xs font-medium text-fg">Condition applicability by model</summary>
              {Object.entries(responsiveness).map(([head, axesOf]) => {
                const inert = Object.entries(axesOf)
                  .filter(([, st]) => st === "fixed")
                  .map(([ax]) => ax);
                const predicts = Object.entries(axesOf)
                  .filter(([, st]) => st === "predicts")
                  .map(([ax]) => ax);
                return (
                  <p key={head} className="text-fg-subtle">
                    <span className="font-mono text-fg">{head}</span>{" "}
                    {inert.length === Object.keys(axesOf).length ? (
                      <span className="text-amber-700">
                        has no condition variation supported by this record;
                        interpret its output within the recorded reference conditions
                      </span>
                    ) : (
                      <>
                        training conditions varied for{" "}
                        {Object.entries(axesOf)
                          .filter(([, st]) => st === "varied")
                          .map(([ax]) => ax)
                          .join(", ") || "no axis"}
                        {inert.length ? `; fixed or unused: ${inert.join(", ")}` : ""}
                        {Object.entries(axesOf).some(([, st]) => st === "unknown")
                          ? `; unknown: ${Object.entries(axesOf).filter(([, st]) => st === "unknown").map(([ax]) => ax).join(", ")}`
                          : ""}
                        {predicts.length ? `; predicts ${predicts.join(", ")}` : ""}
                      </>
                    )}
                  </p>
                );
              })}
              {scan.run_record.condition_note ? (
                <p className="mt-1 text-fg-subtle">{scan.run_record.condition_note}</p>
              ) : null}
            </details>
          ) : null}

          <Quadrants counts={scan.quadrant_counts} />
          <p className="mt-2 text-xs text-fg-muted" data-testid="variant-coverage">
            {scan.variants.length} variants · {scan.variants.filter((v) => v.regulatory !== null).length} with regulatory scores · {scan.variants.filter((v) => v.structural.delta != null).length} with structural deltas · {scan.variants.filter((v) => v.structural.basis === "motif_lost").length} motif-destroying (ranked structurally high without a delta).
            Unclassified variants remain in the table; a missing value is not zero.
          </p>

          {/* Preserve motif events, including those without numerical deltas. */}
          {candidateTotal ? (
            <details
              className="mt-3 rounded-[var(--radius-md)] border border-accent/30 bg-accent/5 p-2"
              data-testid="variant-candidates"
              open
            >
              <summary className="cursor-pointer text-xs font-medium text-fg">
                {candidateTotal} structural candidate
                {candidateTotal === 1 ? "" : "s"} for the selected structural head
              </summary>
              <p className="mt-1 text-[0.65rem] leading-relaxed text-fg-subtle">
                {scan.candidate_note}
              </p>
              <ul className="mt-2 space-y-1 text-[0.7rem]">
                {Object.entries(scan.structural_candidates).flatMap(([state, entries]) =>
                  entries.map((e) => (
                    <li key={`${state}-${e.label}-${e.head}-${e.strand}`} className="flex flex-wrap gap-2">
                      <span className="font-mono text-fg">{e.label}</span>
                      {e.head ? <span className="text-fg-subtle">{e.head}</span> : null}
                      <Badge variant="outline" className="text-[9px]">
                        {state.replace(/_/g, " ")}
                      </Badge>
                      {e.strand ? (
                        <span className="text-fg-subtle">strand {e.strand}</span>
                      ) : null}
                      {e.has_delta && e.delta != null ? (
                        <span className="text-fg-muted">predicted Δ {e.delta.toFixed(3)} {e.units}</span>
                      ) : (
                        <span className="text-fg-subtle">no applicable numerical comparison</span>
                      )}
                      {e.why_listed ? (
                        <span className="basis-full text-fg-subtle">{e.why_listed}</span>
                      ) : null}
                      <span className="text-fg-muted">
                        {e.regulatory_score !== null
                          ? `regulatory ${e.regulatory_score.toFixed(3)}`
                          : "no regulatory record"}
                      </span>
                    </li>
                  )),
                )}
              </ul>
            </details>
          ) : null}

          <p className="mt-2 text-[0.65rem] leading-relaxed text-fg-subtle">
            {scan.quadrant_note} Thresholds: {scan.run_record.threshold_basis}.
          </p>

          {top.length ? (
            <div className="mt-3 min-w-0 max-w-full overflow-x-auto">
              <table className="w-full min-w-[38rem] text-left text-xs">
                <caption className="sr-only">
                  Variants with ranked rows first; unclassified rows retained
                </caption>
                <thead className="text-fg-subtle">
                  <tr>
                    <th scope="col" className="py-1 pr-3 font-normal">Variant</th>
                    <th scope="col" className="py-1 pr-3 font-normal">
                      Δ{units ? ` (${units})` : ""}
                    </th>
                    <th scope="col" className="py-1 pr-3 font-normal">Regulatory</th>
                    <th scope="col" className="py-1 pr-3 font-normal">Rank</th>
                    <th scope="col" className="py-1 font-normal">Quadrant</th>
                  </tr>
                </thead>
                <tbody className="font-mono">
                  {top.map((row) => (
                    <tr key={row.label} className="border-t border-border/50" data-testid="variant-row">
                      <td className="py-1 pr-3 text-fg">
                        {row.label}
                        <span className="ml-1 text-fg-subtle">
                          {row.mutation}
                          {row.structural.strand ? ` · ${row.structural.strand}` : ""}
                        </span>
                      </td>
                      <td className="py-1 pr-3 text-fg">
                        {row.structural.delta != null
                          ? `${row.structural.delta > 0 ? "+" : ""}${row.structural.delta.toFixed(2)}`
                          : row.structural.basis === "motif_lost"
                            ? "motif lost"
                            : "—"}
                      </td>
                      <td className="py-1 pr-3 text-fg">
                        {row.regulatory
                          ? `${row.regulatory.score.toFixed(3)}`
                          : "—"}
                        {row.regulatory ? (
                          <span className="ml-1 text-fg-subtle">{row.regulatory.scorer}</span>
                        ) : null}
                      </td>
                      <td className="py-1 pr-3 text-fg-muted">
                        {row.combined_rank != null ? row.combined_rank.toFixed(2) : "—"}
                      </td>
                      <td className="py-1" title={row.quadrant_reason ?? undefined}>
                        <Badge variant="outline" className="text-[9px]">
                          {QUADRANT_LABEL[row.quadrant]}
                        </Badge>
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          ) : (
            <p className="mt-3 text-xs text-fg-muted">
              No variant rows were returned for this window.
            </p>
          )}
          {scan.variants.length > 25 ? (
            <Button size="sm" variant="ghost" onClick={() => setShowAll(!showAll)}>
              {showAll ? "Show first 25 variants" : `Show all ${scan.variants.length} variants`}
            </Button>
          ) : null}

          <p className="mt-3 flex gap-2 text-[0.7rem] leading-relaxed text-fg-subtle">
            <AlertTriangle className="mt-0.5 h-3 w-3 shrink-0" aria-hidden />
            <span>
              {scan.interpretation_note} {scan.rank_note}
            </span>
          </p>
        </>
      ) : (
        <p className="mt-3 text-xs text-fg-muted">
          Joins the mutation scan of the sequence above to the regulatory score for
          each of the same substitutions. Needs a genomic coordinate, because the
          join key is the variant and not the offset in a pasted string.
        </p>
      )}
    </section>
  );
}
