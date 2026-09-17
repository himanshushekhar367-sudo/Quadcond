/**
 * How does this sequence respond across a buffer range?
 *
 * The interface answered this one point at a time: move a slider, watch a
 * number. Three things that does badly. It gives no view of where the supported
 * region ends. It makes a comparison the reader has to hold in their head. And
 * it puts eight controls in front of twelve heads, of which four respond to any
 * of them — so a slider that moves while a number sits still teaches the reader
 * that the head is *insensitive* to that variable, which is a scientific claim
 * nobody made.
 *
 * This panel replaces the watching with a series, and answers the second
 * question explicitly rather than by omission:
 *
 * - Heads whose training data never varied the axis are **listed as inert with
 *   the reason**, not dropped. A head missing from a chart reads as "not
 *   applicable"; a head listed as inert says which of two very different things
 *   applies — an absence of evidence, or a refusal for this sequence.
 * - Every series carries the **measured** response, not the licensed one.
 *   "The axis varied in training" permits a dependence; it does not show one.
 *   The spread is compared with the head's own error scale, because a chart
 *   auto-scales and a 0.2 °C wobble fills the panel exactly like a 20 °C swing.
 * - Out-of-domain and refused points are drawn as gaps, not interpolated
 *   through. Until v0.4.8 the segmentation checked only whether `y` was null,
 *   so an out-of-domain point that still carried a number was joined into the
 *   path and given an amber dot — the line said "supported response" and the
 *   dot said otherwise. Extrapolation is now an explicit toggle, off by
 *   default, and the response summary the service computes covers in-domain
 *   points only.
 * - **The scan is pinned to what it is a scan of.** `MutationExplorer` has
 *   carried a staleness guard since v0.4.5; this panel had none, so changing
 *   the sequence or the buffer left a completed sweep on screen under a new
 *   heading, and a slow request could land after the user had moved on and
 *   repopulate the panel with the old one. Both are the same failure the
 *   frozen-baseline work was about: a table and its stated reference drifting
 *   apart with nothing saying so.
 */
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { Activity, Download, Loader2 } from "lucide-react";
import { useNA } from "@/lib/na/store";
import {
  QuadcondUnavailable,
  scanConditions,
  type ConditionScan,
  type ConditionSeries,
} from "@/lib/quadcond/client";
import { conditionScanToCsv, download, exportName, toJson } from "@/lib/quadcond/export";
import { Button } from "@/components/ui/button";

/** Axes worth offering, with a default span drawn from what gets measured. */
const AXES: { id: string; label: string; min: number; max: number; unit: string }[] = [
  { id: "k", label: "K⁺", min: 0, max: 150, unit: "mM" },
  { id: "na", label: "Na⁺", min: 0, max: 150, unit: "mM" },
  { id: "mg", label: "Mg²⁺", min: 0, max: 10, unit: "mM" },
  { id: "ph", label: "pH", min: 4.5, max: 8.0, unit: "" },
];

const VERDICT_TONE: Record<string, string> = {
  responds: "text-emerald-700",
  // `flat` is retired. It read as a measured insensitivity to the axis, which
  // is a claim nobody made: what was computed is a response range beside a
  // marginal residual half-width.
  below_error_scale: "text-amber-700",
  unknown_scale: "text-fg-muted",
  refused_throughout: "text-amber-700",
  insufficient_points: "text-fg-muted",
  insufficient_supported_points: "text-amber-700",
};

const VERDICT_LABEL: Record<string, string> = {
  responds: "responds",
  below_error_scale: "below this head's error scale",
  unknown_scale: "no error scale to judge against",
  refused_throughout: "refused at every point",
  insufficient_points: "too few points",
  insufficient_supported_points: "too few supported points",
};

/**
 * A minimal inline plot.
 *
 * Deliberately not a charting library: the only thing being drawn is one
 * series against one axis, and the interesting design decisions are about what
 * *not* to draw — points the head refused are gaps, and the y-scale is annotated
 * with the head's error scale so the eye cannot mistake noise for a response.
 */
function Sparkline({
  series,
  showUnsupported,
}: {
  series: ConditionSeries;
  showUnsupported: boolean;
}) {
  const key = series.observed_response.quantity ?? "value";
  const pts = series.points.map((p) => ({
    x: p.value_of_axis,
    y: (key === "value" ? p.value : p.probability) ?? null,
    ok: p.applicability.state === "in_domain",
    refused: p.applicability.state === "refused",
  }));
  // The supported points set the scale. Letting an extrapolated tail stretch
  // the y-axis compresses the part of the curve the head can actually answer
  // for into a flat-looking band.
  const ys = pts
    .filter((p) => p.ok)
    .map((p) => p.y)
    .filter((y): y is number => y !== null);
  if (ys.length < 2) return null;

  const xs = pts.map((p) => p.x);
  const [x0, x1] = [Math.min(...xs), Math.max(...xs)];
  const [y0, y1] = [Math.min(...ys), Math.max(...ys)];
  const pad = (y1 - y0) * 0.15 || 1;
  const W = 260;
  const H = 64;
  const sx = (x: number) => ((x - x0) / (x1 - x0 || 1)) * (W - 8) + 4;
  const sy = (y: number) => H - 6 - ((y - (y0 - pad)) / (y1 + pad - (y0 - pad) || 1)) * (H - 12);

  // Segments, not one path. A point is a break in the line when there is no
  // value **or when the head reports the query out of its domain** — the second
  // half was missing, so an extrapolated point with a number in it was joined
  // into the path and merely coloured amber. A drawn line asserts a supported
  // response between its endpoints; a dot colour does not retract that.
  const drawable = (p: { y: number | null; ok: boolean; refused: boolean }) =>
    !p.refused && p.y !== null && Number.isFinite(p.y) && (p.ok || showUnsupported);
  const segments: string[] = [];
  let run: string[] = [];
  for (const p of pts) {
    if (!drawable(p)) {
      if (run.length > 1) segments.push(run.join(" "));
      run = [];
      continue;
    }
    run.push(`${run.length ? "L" : "M"}${sx(p.x).toFixed(1)},${sy(p.y!).toFixed(1)}`);
  }
  if (run.length > 1) segments.push(run.join(" "));

  // The head's own typical error, drawn to the same scale as the data. If the
  // whole series fits inside this bar, its range is below the marginal error
  // scale. This does not test the uncertainty of a difference.
  const scale = series.observed_response.error_scale ?? null;
  const errBar = scale !== null ? (scale / (y1 + pad - (y0 - pad) || 1)) * (H - 12) : null;

  return (
    <svg
      viewBox={`0 0 ${W} ${H}`}
      className="mt-2 h-16 w-full"
      role="img"
      aria-label={`${series.head} across the scanned range: ${series.observed_response.note ?? ""}`}
    >
      {errBar !== null ? (
        <g>
          <rect
            x={W - 14}
            y={H / 2 - errBar / 2}
            width={4}
            height={Math.max(errBar, 1)}
            className="fill-fg-subtle/40"
          />
          <title>Typical error for this head, to scale</title>
        </g>
      ) : null}
      {segments.map((d, i) => (
        <path key={i} d={d} className="fill-none stroke-accent" strokeWidth={1.5} />
      ))}
      {pts.map((p, i) =>
        !drawable(p) ? null : (
          <circle
            key={i}
            cx={sx(p.x)}
            cy={sy(p.y!)}
            data-testid="condition-dot"
            data-supported={p.ok}
            r={2}
            className={p.ok ? "fill-accent" : "fill-amber-400"}
          >
            <title>
              {p.ok
                ? `${p.x}: ${p.y}`
                : `${p.x}: ${p.y} — outside this head's applicability domain, ` +
                  `excluded from the response summary`}
            </title>
          </circle>
        ),
      )}
    </svg>
  );
}

export function ConditionResponse() {
  const sequence = useNA((s) => s.sequence);
  const conditions = useNA((s) => s.conditions);
  const backend = useNA((s) => s.backend);

  const [axis, setAxis] = useState("k");
  const [scan, setScan] = useState<ConditionScan | null>(null);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [showUnsupported, setShowUnsupported] = useState(false);

  const spec = useMemo(() => AXES.find((a) => a.id === axis)!, [axis]);

  /**
   * What this sweep is a sweep *of*.
   *
   * Same shape as `MutationExplorer`'s guard, and here for the same reason: it
   * remembers the inputs that were sent rather than re-deriving them, because a
   * guard that re-derives what it is checking is a second implementation of the
   * thing it is guarding and will drift. (That drift is not hypothetical — the
   * first version of the mutation guard reimplemented the backend's sequence
   * cleaning, disagreed with it on one character, and discarded every scan in
   * the tick it arrived.)
   */
  const conditionKey = JSON.stringify(conditions);
  const [pin, setPin] = useState<
    { sequence: string; conditions: string; axis: string } | null
  >(null);
  const stale =
    pin !== null &&
    (pin.sequence !== sequence || pin.conditions !== conditionKey || pin.axis !== axis);

  useEffect(() => {
    // A sweep under a different sequence, buffer or axis is not slightly out of
    // date — it is about a different question. Dropped rather than left on
    // screen under a new heading.
    if (stale) {
      setScan(null);
      setPin(null);
      setError(null);
    }
  }, [stale]);

  /**
   * Which request the panel is currently waiting for.
   *
   * A sweep is several seconds of real model work. Without this, changing the
   * axis and sweeping again could land the *first* response last and repopulate
   * the panel with a scan of something the user had already moved off — a
   * result that is not merely stale but silently mislabelled by the header
   * around it.
   */
  const runId = useRef(0);

  const run = useCallback(async () => {
    const id = ++runId.current;
    const pinned = { sequence, conditions: conditionKey, axis };
    setBusy(true);
    setError(null);
    try {
      const result = await scanConditions(
        sequence,
        axis,
        { min: spec.min, max: spec.max, points: 12 },
        conditions,
      );
      if (id !== runId.current) return;      // a newer request has been issued
      setScan(result);
      setPin(pinned);
    } catch (err) {
      if (id !== runId.current) return;
      setError(err instanceof QuadcondUnavailable ? err.message : String(err));
      setScan(null);
    } finally {
      if (id === runId.current) setBusy(false);
    }
  }, [sequence, axis, spec, conditions, conditionKey]);

  if (backend !== "connected") return null;

  return (
    <section
      className="rounded-[var(--radius-lg)] border border-border bg-surface p-4"
      data-testid="condition-response"
    >
      <header className="flex flex-wrap items-center justify-between gap-2">
        <div>
          <h3 className="flex items-center gap-2 text-sm font-medium text-fg">
            <Activity className="h-4 w-4" aria-hidden />
            Condition response
          </h3>
          <p className="text-xs text-fg-subtle">
            One axis, swept — with the heads that cannot answer named rather than
            hidden.
          </p>
        </div>
        {scan ? (
          <div className="flex gap-1">
            <Button
              size="sm"
              variant="ghost"
              onClick={() =>
                download(
                  exportName(`condition-${scan.axis}`, scan.sequence, "csv"),
                  conditionScanToCsv(scan),
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
                  exportName(`condition-${scan.axis}`, scan.sequence, "json"),
                  toJson(scan),
                  "application/json",
                )
              }
            >
              JSON
            </Button>
          </div>
        ) : null}
      </header>

      <div className="mt-3 flex flex-wrap items-center gap-1.5">
        {AXES.map((a) => (
          <button
            key={a.id}
            type="button"
            onClick={() => {
              setAxis(a.id);
              setScan(null);
              setPin(null);
            }}
            data-testid="condition-axis"
            className={`rounded-full border px-2.5 py-1 text-xs transition ${
              axis === a.id
                ? "border-accent bg-accent/15 text-fg"
                : "border-border text-fg-muted hover:text-fg"
            }`}
          >
            {a.label}
          </button>
        ))}
        <Button size="sm" onClick={run} disabled={busy} className="ml-auto">
          {busy ? <Loader2 className="h-3.5 w-3.5 animate-spin" aria-hidden /> : null}
          Sweep {spec.min}–{spec.max} {spec.unit}
        </Button>
      </div>

      {error ? (
        <p className="mt-3 rounded-[var(--radius-md)] border border-amber-500/40 bg-amber-500/10 p-2 text-xs text-amber-700">
          {error}
        </p>
      ) : null}

      {scan ? (
        <>
          {/* The frozen context, read off the scan rather than off live state.
              These are the lines that could previously disagree with the chart
              above them. */}
          <dl
            className="mt-3 grid grid-cols-2 gap-x-4 gap-y-1 rounded-[var(--radius-md)] border border-border bg-bg/40 p-3 text-xs sm:grid-cols-4"
            data-testid="condition-frozen-context"
          >
            <div className="col-span-2">
              <dt className="text-fg-subtle">Swept sequence</dt>
              <dd className="font-mono break-all text-fg">{scan.sequence}</dd>
            </div>
            <div>
              <dt className="text-fg-subtle">Buffer</dt>
              <dd className="text-fg">
                {Object.entries(scan.base_condition)
                  .map(([k, v]) => `${k} ${v}`)
                  .join(" · ")}
              </dd>
            </div>
            <div>
              <dt className="text-fg-subtle">Model</dt>
              <dd className="font-mono text-fg">{scan.model_version}</dd>
            </div>
            <div className="col-span-2 sm:col-span-4">
              <dt className="text-fg-subtle">Reverse complement (strand −)</dt>
              <dd className="font-mono break-all text-fg-muted">{scan.complement}</dd>
            </div>
            {scan.condition_imputed_fields?.length ? (
              <div className="col-span-2 sm:col-span-4">
                <dt className="text-fg-subtle">Defaulted, not specified</dt>
                <dd className="text-fg-muted">
                  {scan.condition_imputed_fields.join(", ")}
                </dd>
              </div>
            ) : null}
          </dl>

          <label className="mt-2 flex items-center gap-2 text-[0.7rem] text-fg-muted">
            <input
              type="checkbox"
              checked={showUnsupported}
              onChange={(e) => setShowUnsupported(e.target.checked)}
              data-testid="condition-show-unsupported"
            />
            Draw extrapolated points. Off by default: a drawn line asserts a
            supported response between its endpoints, and these points are
            outside the head&rsquo;s domain and excluded from every summary
            below.
          </label>

          <div className="mt-3 grid gap-3 sm:grid-cols-2">
            {Object.entries(scan.series).map(([name, series]) => (
              <div
                key={name}
                className="rounded-[var(--radius-md)] border border-border bg-bg/40 p-3"
                data-testid="condition-series"
              >
                <div className="flex items-baseline justify-between gap-2">
                  <span className="font-mono text-xs text-fg">
                    {name}
                    {/* Which strand this head was asked about. Two cards in one
                        grid are otherwise about two different molecules. */}
                    <span className="ml-1 text-fg-subtle">strand {series.strand}</span>
                  </span>
                  <span
                    className={`text-[0.65rem] ${VERDICT_TONE[series.observed_response.verdict] ?? ""}`}
                  >
                    {VERDICT_LABEL[series.observed_response.verdict] ??
                      series.observed_response.verdict.replace(/_/g, " ")}
                  </span>
                </div>
                <Sparkline series={series} showUnsupported={showUnsupported} />
                <p className="mt-1 text-[0.65rem] leading-relaxed text-fg-subtle">
                  {series.observed_response.note ??
                    "No response summary for this head."}
                </p>
                {series.observed_response.n_out_of_domain ? (
                  <p className="mt-1 text-[0.65rem] text-amber-700/90">
                    {series.observed_response.n_out_of_domain} of{" "}
                    {series.observed_response.n_points} points are outside this
                    head&rsquo;s domain and are not in the summary above.
                  </p>
                ) : null}
              </div>
            ))}
          </div>

          {/* The half of the answer a chart cannot show. */}
          {scan.inert_heads.length ? (
            <details className="mt-3 rounded-[var(--radius-md)] border border-border p-2">
              <summary className="cursor-pointer text-xs text-fg-muted">
                {scan.inert_heads.length} head
                {scan.inert_heads.length === 1 ? "" : "s"} cannot respond to{" "}
                {scan.axis_label} — and why
              </summary>
              <ul className="mt-2 space-y-1.5 text-[0.7rem]">
                {scan.inert_heads.map((h) => (
                  <li key={h.head}>
                    <span className="font-mono text-fg">{h.head}</span>{" "}
                    <span className="text-fg-subtle">{h.reason}</span>
                  </li>
                ))}
              </ul>
            </details>
          ) : null}

          {scan.summary_note ? (
            <p className="mt-2 text-[0.65rem] text-fg-subtle">{scan.summary_note}</p>
          ) : null}

          {scan.measured_neighbours.length ? (
            <p className="mt-2 text-[0.65rem] text-fg-subtle">
              {scan.neighbour_note}
            </p>
          ) : null}
        </>
      ) : null}
    </section>
  );
}
