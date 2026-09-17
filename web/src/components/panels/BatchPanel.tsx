/**
 * Many sequences at once, and a file you can defend six months later.
 *
 * The interface had one sequence box. That is fine for looking at a structure
 * and useless for the thing this tool is actually for — deciding which of forty
 * candidates to put on an instrument. Screenshots do not survive review, and a
 * table of melting temperatures with no head, no model version and no buffer is
 * a table of plausible numbers.
 *
 * Two things this panel takes seriously:
 *
 * **Identifiers survive.** FASTA headers come back verbatim on every row. A
 * result that renumbers the caller's sequences has to be re-matched by hand, and
 * hand-matching is where a row acquires the wrong label.
 *
 * **Exclusions are shown.** A batch of 200 that returns 197 rows with no note is
 * a batch whose three missing sequences are noticed after the figure is drawn,
 * if at all. Anything the service declined is listed with its reason, in the
 * panel and in the export.
 *
 * **And every supplied record is reconciled.** Until v0.4.8 the panel's label
 * promised "FASTA or one sequence per line" and the parser delivered neither
 * reliably: a header with no sequence produced no record at all rather than an
 * excluded one, bare lines were concatenated into a single sequence, and the
 * scored count was taken over identifiers, so two records sharing a name were
 * reported as one. Three inputs could become two results and an empty exclusion
 * list. The run record now counts by input index and states that every record
 * landed in exactly one of results or exclusions.
 */
import { useCallback, useMemo, useState } from "react";
import { Download, Loader2, Rows3 } from "lucide-react";
import { useNA } from "@/lib/na/store";
import { batch, QuadcondUnavailable, type BatchResult } from "@/lib/quadcond/client";
import { batchToCsv, download, exportName, toJson } from "@/lib/quadcond/export";
import { Button } from "@/components/ui/button";

const PLACEHOLDER = `>tel22
AGGGTTAGGGTTAGGGTTAGGG
>myc-pu27
TGGGGAGGGTGGGGAGGGTGGGGAAGG
>ckit1
AGGGAGGGCGCTGGGAGGAGGG`;

/** Which heads a row actually carries, so the table has columns to draw. */
function headsIn(result: BatchResult): string[] {
  const seen = new Set<string>();
  for (const row of result.results) {
    for (const name of Object.keys(
      (row.predictions ?? {}) as Record<string, unknown>,
    )) {
      seen.add(name);
    }
  }
  return [...seen].sort();
}

export function BatchPanel() {
  const conditions = useNA((s) => s.conditions);
  const backend = useNA((s) => s.backend);

  const [text, setText] = useState("");
  const [result, setResult] = useState<BatchResult | null>(null);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const heads = useMemo(() => (result ? headsIn(result) : []), [result]);

  const run = useCallback(async () => {
    setBusy(true);
    setError(null);
    try {
      // `auto`: FASTA when any line starts with ">", one sequence per line
      // otherwise — which is what the label under this box promises.
      setResult(
        await batch({ fasta: text || PLACEHOLDER, parse_mode: "auto" }, conditions),
      );
    } catch (err) {
      setError(err instanceof QuadcondUnavailable ? err.message : String(err));
      setResult(null);
    } finally {
      setBusy(false);
    }
  }, [text, conditions]);

  if (backend !== "connected") return null;

  return (
    <section
      className="rounded-[var(--radius-lg)] border border-border bg-surface p-4"
      data-testid="batch-panel"
    >
      <header className="flex flex-wrap items-center justify-between gap-2">
        <div>
          <h3 className="flex items-center gap-2 text-sm font-medium text-fg">
            <Rows3 className="h-4 w-4" aria-hidden />
            Batch
          </h3>
          <p className="text-xs text-fg-subtle">
            FASTA, or one sequence per line. Every record comes back as a row or
            as a stated exclusion.
          </p>
        </div>
        <div className="flex gap-1">
          <Button size="sm" onClick={run} disabled={busy}>
            {busy ? <Loader2 className="h-3.5 w-3.5 animate-spin" aria-hidden /> : null}
            Run
          </Button>
          {result ? (
            <>
              <Button
                size="sm"
                variant="ghost"
                onClick={() =>
                  download(exportName("batch", "run", "csv"), batchToCsv(result, heads))
                }
              >
                <Download className="h-3.5 w-3.5" aria-hidden /> CSV
              </Button>
              <Button
                size="sm"
                variant="ghost"
                onClick={() =>
                  download(
                    exportName("batch", "run", "json"),
                    toJson(result),
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

      <label className="mt-3 block">
        <span className="sr-only">
          FASTA, or one sequence per line. A header with no sequence is reported
          as an exclusion, not dropped.
        </span>
        <textarea
          value={text}
          onChange={(e) => setText(e.target.value)}
          placeholder={PLACEHOLDER}
          rows={6}
          spellCheck={false}
          data-testid="batch-input"
          className="w-full rounded-[var(--radius-md)] border border-border bg-bg/60 p-2 font-mono text-xs text-fg placeholder:text-fg-subtle focus:border-accent focus:outline-none"
        />
      </label>

      {error ? (
        <p className="mt-2 rounded-[var(--radius-md)] border border-amber-500/40 bg-amber-500/10 p-2 text-xs text-amber-700">
          {error}
        </p>
      ) : null}

      {result ? (
        <>
          <div className="mt-3 min-w-0 max-w-full overflow-x-auto">
            <table className="w-full min-w-[30rem] text-left text-xs">
              <caption className="sr-only">Batch predictions by sequence and head</caption>
              <thead className="text-fg-subtle">
                <tr>
                  <th scope="col" className="py-1 pr-3 font-normal">ID</th>
                  {heads.map((h) => (
                    <th key={h} scope="col" className="py-1 pr-3 font-normal">
                      {h}
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody className="font-mono">
                {result.results.map((row, i) => {
                  const preds = (row.predictions ?? {}) as Record<
                    string,
                    Record<string, unknown>
                  >;
                  return (
                    <tr key={`${row.id}-${i}`} className="border-t border-border/50">
                      <td className="py-1 pr-3 text-fg">{row.id}</td>
                      {heads.map((h) => {
                        const e = preds[h];
                        if (!e) return <td key={h} className="py-1 pr-3 text-fg-subtle">—</td>;
                        // A refused head has no value key. Rendering "refused"
                        // rather than an empty cell keeps it distinguishable
                        // from a head that simply was not requested.
                        if (e.refused) {
                          return (
                            <td key={h} className="py-1 pr-3 text-amber-700" title={String(e.refusal_reason ?? "")}>
                              refused
                            </td>
                          );
                        }
                        const ap = (e.applicability ?? {}) as Record<string, unknown>;
                        const v = e.value ?? e.probability ?? e.argmax;
                        return (
                          <td
                            key={h}
                            className={`py-1 pr-3 ${ap.in_domain === false ? "text-amber-700" : "text-fg"}`}
                            title={
                              ap.in_domain === false
                                ? `out of domain: ${((ap.warnings as string[]) ?? []).join(" | ")}`
                                : undefined
                            }
                          >
                            {typeof v === "number" ? v.toFixed(2) : String(v ?? "—")}
                          </td>
                        );
                      })}
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </div>

          {result.excluded.length ? (
            <div className="mt-3 rounded-[var(--radius-md)] border border-amber-500/40 bg-amber-500/10 p-2">
              <p className="text-xs font-medium text-amber-700">
                {result.excluded.length} input
                {result.excluded.length === 1 ? "" : "s"} not scored
              </p>
              <ul className="mt-1 space-y-0.5 text-[0.7rem] text-amber-700/90">
                {result.excluded.map((e) => (
                  <li key={`${e.record_index ?? e.id}`}>
                    <span className="font-mono">{e.id}</span>
                    {e.record_index !== undefined ? (
                      <span className="text-amber-700/60"> (record {e.record_index + 1})</span>
                    ) : null}{" "}
                    — {e.reason}
                  </li>
                ))}
              </ul>
            </div>
          ) : null}

          {/* The count line, stated rather than left to be inferred from a row
              tally. "3 in, 2 scored, 1 excluded" is checkable at a glance; three
              inputs silently becoming two rows is not. */}
          <p className="mt-3 text-xs text-fg-muted" data-testid="batch-reconciliation">
            {result.run_record.n_sequences_in} record
            {result.run_record.n_sequences_in === 1 ? "" : "s"} in ·{" "}
            {result.run_record.n_sequences_scored} scored ·{" "}
            {result.run_record.n_sequences_excluded ?? result.excluded.length} excluded ·{" "}
            {result.run_record.n_rows} rows
            {result.run_record.duplicate_identifiers?.length ? (
              <span className="text-amber-700">
                {" "}
                · {result.run_record.duplicate_identifiers.length} identifier
                {result.run_record.duplicate_identifiers.length === 1 ? "" : "s"} used
                more than once ({result.run_record.duplicate_identifiers.join(", ")}) —
                counted as separate sequences
              </span>
            ) : null}
          </p>

          <details className="mt-3 rounded-[var(--radius-md)] border border-border p-2">
            <summary className="cursor-pointer text-xs text-fg-muted">
              Run record — what the export carries
            </summary>
            <dl className="mt-2 grid grid-cols-2 gap-x-3 gap-y-1 text-[0.7rem]">
              {Object.entries(result.run_record).map(([k, v]) => (
                <div key={k} className="contents">
                  <dt className="text-fg-subtle">{k}</dt>
                  <dd className="truncate font-mono text-fg" title={JSON.stringify(v)}>
                    {typeof v === "object" ? JSON.stringify(v) : String(v)}
                  </dd>
                </div>
              ))}
            </dl>
          </details>
        </>
      ) : null}
    </section>
  );
}
