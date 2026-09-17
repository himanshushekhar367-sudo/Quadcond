/**
 * Which single substitution is worth measuring next?
 *
 * A branch of this project shipped a panel with this name that answered the
 * question badly enough to be worse than not answering it. It called the local
 * heuristic energy function — the one that cannot distinguish K⁺ from Na⁺ and
 * scores R² 0.074 against measured folding free energy — ranked substitutions
 * by its ΔΔG, and labelled anything past a fixed magnitude "significant". It
 * also recomputed the displayed wild-type baseline when a mutant was selected
 * while keeping the previous scan's rows on screen, so the table and the
 * reference it was supposedly relative to could disagree with nothing saying so.
 *
 * This version asks QuadCond, and the differences from that one are all the
 * same difference: a number here has to be about something.
 *
 * - **The head that owns the quantity answers.** ΔTm comes from `g4_tm`,
 *   ΔpH_T from `im_pht`. There is no cross-structure energy difference, because
 *   a difference between two structures' energies is not a property of either.
 * - **Each head is asked about its own strand.** A G4 and an i-motif at one
 *   duplex position sit on opposite strands; one edit is propagated to both.
 * - **The baseline is frozen** into the scan and displayed from it, so the
 *   table and its stated reference cannot drift apart.
 * - **Changing conditions invalidates the scan** rather than leaving stale rows
 *   under a new buffer.
 * - **Motif loss and out-of-domain are states, not blanks.** A substitution
 *   that destroys the motif produces no delta at all — no key, no zero.
 * - **Nothing is called significant.** Rows are ranked by magnitude and carry
 *   the paired spread; there is no null distribution here to test against.
 */
import { useCallback, useEffect, useMemo, useState } from "react";
import { AlertTriangle, Download, FlaskConical, Loader2 } from "lucide-react";
import { useNA } from "@/lib/na/store";
import {
  QuadcondUnavailable,
  scanMutations,
  type MutationCell,
  type MutationScan,
} from "@/lib/quadcond/client";
import { download, exportName, mutationScanToCsv, toJson } from "@/lib/quadcond/export";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";

const MOTIF_LABEL: Record<MutationCell["motif"]["state"], string> = {
  motif_retained: "motif kept",
  motif_lost: "motif destroyed",
  motif_gained: "motif created",
  motif_count_changed: "motif count changed",
  no_motif: "no motif on this strand",
};

/** Rows with a number first, then the ones whose finding is a state. */
function rank(scan: MutationScan, head: string) {
  const withDelta = scan.substitutions.filter((r) => r.heads[head]?.delta !== undefined);
  const without = scan.substitutions.filter((r) => r.heads[head]?.delta === undefined);
  withDelta.sort(
    (a, b) => Math.abs(b.heads[head]!.delta!) - Math.abs(a.heads[head]!.delta!),
  );
  return { withDelta, without };
}

export function MutationExplorer() {
  const sequence = useNA((s) => s.sequence);
  const conditions = useNA((s) => s.conditions);
  const backend = useNA((s) => s.backend);

  const [scan, setScan] = useState<MutationScan | null>(null);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [head, setHead] = useState<string | null>(null);

  /**
   * What this scan is a scan *of*.
   *
   * The staleness guard compares the inputs that were sent with the inputs as
   * they are now — not a re-derivation of them. The first version of this
   * reimplemented the backend's sequence cleaning in the client to compare
   * against `scan.wild_type.sequence`, and got it subtly wrong: Python masks an
   * unrecognised character to `N` and keeps the length, the client regex
   * deleted it. Any sequence with a stray character produced two different
   * strings, the guard read that as "the user has moved on", and **every scan
   * was discarded in the same tick it arrived** — no rows, no error, the panel
   * back at its start state as though the button had never been pressed.
   *
   * A guard that re-derives what it is checking is a second implementation of
   * the thing it is guarding, and it will drift. This one remembers.
   */
  const [pin, setPin] = useState<{ sequence: string; conditions: string } | null>(null);
  const conditionKey = JSON.stringify(conditions);
  const stale = pin !== null && (pin.sequence !== sequence || pin.conditions !== conditionKey);

  useEffect(() => {
    // A result under a different sequence or buffer is not slightly out of
    // date, it is about a different question. Dropped rather than left on
    // screen under a new heading — which is what let the previous panel show a
    // table and a baseline that disagreed.
    if (stale) {
      setScan(null);
      setPin(null);
      setError(null);
    }
  }, [stale]);

  const run = useCallback(async () => {
    setBusy(true);
    setError(null);
    try {
      const result = await scanMutations(sequence, conditions);
      setScan(result);
      setPin({ sequence, conditions: conditionKey });
      setHead((h) => (h && result.heads[h] ? h : (Object.keys(result.heads)[0] ?? null)));
    } catch (err) {
      setError(
        err instanceof QuadcondUnavailable
          ? err.message
          : `scan failed: ${String(err)}`,
      );
      setScan(null);
    } finally {
      setBusy(false);
    }
  }, [sequence, conditions, conditionKey]);

  const ranked = useMemo(
    () => (scan && head ? rank(scan, head) : null),
    [scan, head],
  );

  if (backend !== "connected") {
    return (
      <div className="rounded-[var(--radius-lg)] border border-border bg-surface p-4">
        <p className="text-sm text-fg-muted">
          Mutation scanning needs QuadCond. The local path builds geometry and
          carries no numbers, so there is nothing here to compare — and the
          heuristic that used to answer this offline could not tell potassium
          from sodium.
        </p>
      </div>
    );
  }

  const units = head && scan ? (scan.heads[head]?.units ?? "") : "";

  return (
    <section
      className="rounded-[var(--radius-lg)] border border-border bg-surface p-4"
      data-testid="mutation-explorer"
    >
      <header className="flex flex-wrap items-center justify-between gap-2">
        <div>
          <h3 className="flex items-center gap-2 text-sm font-medium text-fg">
            <FlaskConical className="h-4 w-4" aria-hidden />
            Mutation scan
          </h3>
          <p className="text-xs text-fg-subtle">
            Every single substitution, against a frozen wild type.
          </p>
        </div>
        <div className="flex items-center gap-2">
          <Button size="sm" onClick={run} disabled={busy} data-testid="run-mutation-scan">
            {busy ? <Loader2 className="h-3.5 w-3.5 animate-spin" aria-hidden /> : null}
            {scan ? "Re-run" : "Run scan"}
          </Button>
          {scan ? (
            <>
              <Button
                size="sm"
                variant="ghost"
                onClick={() =>
                  download(
                    exportName("mutation-scan", scan.wild_type.sequence, "csv"),
                    mutationScanToCsv(scan),
                  )
                }
              >
                <Download className="h-3.5 w-3.5" aria-hidden /> CSV
              </Button>
              <Button
                size="sm"
                variant="ghost"
                onClick={() =>
                  download(
                    exportName("mutation-scan", scan.wild_type.sequence, "json"),
                    toJson(scan),
                    "application/json",
                  )
                }
              >
                JSON
              </Button>
            </>
          ) : null}
        </div>
      </header>

      {error ? (
        <p className="mt-3 rounded-[var(--radius-md)] border border-amber-500/40 bg-amber-500/10 p-2 text-xs text-amber-700">
          {error}
        </p>
      ) : null}

      {scan ? (
        <>
          {/* The frozen baseline, read off the scan rather than off live state.
              This is the line that could previously disagree with the table. */}
          <dl className="mt-3 grid grid-cols-2 gap-x-4 gap-y-1 rounded-[var(--radius-md)] border border-border bg-bg/40 p-3 text-xs sm:grid-cols-4">
            <div>
              <dt className="text-fg-subtle">Wild type</dt>
              <dd className="font-mono break-all text-fg">{scan.wild_type.sequence}</dd>
            </div>
            <div>
              <dt className="text-fg-subtle">Buffer</dt>
              <dd className="text-fg">{scan.wild_type.condition}</dd>
            </div>
            <div>
              <dt className="text-fg-subtle">Model</dt>
              <dd className="font-mono text-fg">{scan.model_version}</dd>
            </div>
            <div>
              <dt className="text-fg-subtle">Substitutions</dt>
              <dd className="text-fg">{scan.substitutions.length}</dd>
            </div>
            {scan.wild_type.condition_imputed_fields.length ? (
              <div className="col-span-2 sm:col-span-4">
                <dt className="text-fg-subtle">Defaulted, not specified</dt>
                <dd className="text-fg-muted">
                  {scan.wild_type.condition_imputed_fields.join(", ")}
                </dd>
              </div>
            ) : null}
          </dl>

          <div className="mt-3 flex flex-wrap gap-1.5">
            {Object.entries(scan.heads).map(([name, info]) => (
              <button
                key={name}
                type="button"
                onClick={() => setHead(name)}
                data-testid="mutation-head-tab"
                className={`rounded-full border px-2.5 py-1 text-xs transition ${
                  head === name
                    ? "border-accent bg-accent/15 text-fg"
                    : "border-border text-fg-muted hover:text-fg"
                }`}
              >
                {name}
                <span className="ml-1 text-fg-subtle">
                  strand {scan.wild_type.strand_by_head[name]}
                </span>
              </button>
            ))}
          </div>

          {head && ranked ? (
            <>
              <p className="mt-3 text-xs text-fg-subtle">{scan.ranking_note}</p>

              <div className="mt-2 min-w-0 max-w-full overflow-x-auto">
                <table className="w-full min-w-[34rem] text-left text-xs">
                  <caption className="sr-only">
                    Predicted change per substitution for {head}
                  </caption>
                  <thead className="text-fg-subtle">
                    <tr>
                      <th scope="col" className="py-1 pr-3 font-normal">Mutation</th>
                      <th scope="col" className="py-1 pr-3 font-normal">
                        Δ{units ? ` (${units})` : ""}
                      </th>
                      <th scope="col" className="py-1 pr-3 font-normal">
                        paired SD
                      </th>
                      <th scope="col" className="py-1 pr-3 font-normal">Mutant</th>
                      <th scope="col" className="py-1 font-normal">State</th>
                    </tr>
                  </thead>
                  <tbody className="font-mono">
                    {ranked.withDelta.map((row) => {
                      const c = row.heads[head]!;
                      return (
                        <tr key={row.label} className="border-t border-border/50">
                          <td className="py-1 pr-3 text-fg">{row.label}</td>
                          <td className="py-1 pr-3 text-fg">
                            {c.delta! > 0 ? "+" : ""}
                            {c.delta!.toFixed(2)}
                          </td>
                          <td className="py-1 pr-3 text-fg-muted">
                            ±{c.delta_paired_sd!.toFixed(2)}
                          </td>
                          <td className="py-1 pr-3 text-fg-muted">
                            {c.mutant_value!.toFixed(2)}
                          </td>
                          <td className="py-1">
                            <Badge
                              variant={
                                c.mutant.state === "in_domain" ? "outline" : "default"
                              }
                            >
                              {c.mutant.state.replace(/_/g, " ")}
                            </Badge>
                          </td>
                        </tr>
                      );
                    })}
                  </tbody>
                </table>
              </div>

              {/* Not an appendix of failures. "The motif is gone" is often the
                  strongest result a scan produces, and it has no delta by
                  construction — so it gets its own list rather than a blank
                  cell in the one above. */}
              {ranked.without.length ? (
                <details className="mt-3 rounded-[var(--radius-md)] border border-border p-2">
                  <summary className="cursor-pointer text-xs text-fg-muted">
                    {ranked.without.length} substitution
                    {ranked.without.length === 1 ? "" : "s"} with no comparable
                    number — and why
                  </summary>
                  <ul className="mt-2 space-y-1 text-xs">
                    {ranked.without.map((row) => {
                      const c = row.heads[head]!;
                      return (
                        <li key={row.label} className="flex gap-2">
                          <span className="font-mono text-fg">{row.label}</span>
                          <span className="text-fg-muted">
                            {MOTIF_LABEL[c.motif.state]}
                            {c.mutant.state !== "in_domain"
                              ? ` · ${c.mutant.state.replace(/_/g, " ")}`
                              : ""}
                          </span>
                        </li>
                      );
                    })}
                  </ul>
                </details>
              ) : null}

              <p className="mt-3 flex gap-2 text-[0.7rem] leading-relaxed text-fg-subtle">
                <AlertTriangle className="mt-0.5 h-3 w-3 shrink-0" aria-hidden />
                <span>
                  {/* Read from the top level, where the service actually sends
                      it. This used to read a per-cell `delta_note` the service
                      stopped emitting in v0.4.6 — it found undefined every time
                      and fell back to a one-line paraphrase that dropped the
                      part that matters: the paired spread is not a validated
                      interval on the difference, and mutation-effect prediction
                      has not been validated here as its own task. */}
                  {scan.delta_note}
                </span>
              </p>
            </>
          ) : null}
        </>
      ) : (
        <p className="mt-3 text-xs text-fg-muted">
          Runs {sequence.replace(/[^ACGTUacgtu]/g, "").length * 3} predictions per
          head — every base changed to each of the other three, compared with the
          wild type under the buffer set above.
        </p>
      )}
    </section>
  );
}
