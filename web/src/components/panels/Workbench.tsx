/**
 * The four workspaces, and where the twelve heads actually live.
 *
 * Twelve heads on one screen implies twelve comparable numbers. They are four
 * different kinds of statement: what the sequence supports, how it moves with
 * buffer, what an antibody found at a locus, and what another model said. Put
 * side by side in one list they are read as one list.
 *
 * The grouping is **discovered from `/info`**, not typed in here. A list of head
 * names in a component is a list that goes stale the first time a head is added
 * and has no way to know it has — which is how the interface came to show eleven
 * heads while the model carried twelve. The backend derives the groups from each
 * head's own claim basis (`quadcond/capabilities.py`), so a new head cannot
 * appear without being classified.
 *
 * Empty groups are rendered rather than dropped. A viewer showing three tabs
 * when the backend has four kinds of head has quietly hidden a category of
 * claim.
 */
import { useEffect, useState } from "react";
import { ChevronRight, Info } from "lucide-react";
import { useNA } from "@/lib/na/store";
import { info, type HeadInfo, type ServiceInfo } from "@/lib/quadcond/client";
import { MutationExplorer } from "./MutationExplorer";
import { VariantPanel } from "./VariantPanel";
import { ConditionResponse } from "./ConditionResponse";
import { BatchPanel } from "./BatchPanel";
import { Tabs, TabsList, TabsTrigger, TabsContent } from "@/components/ui/tabs";

const SEMANTICS_TONE: Record<string, string> = {
  biophysical: "text-emerald-700",
  genomic_proxy: "text-amber-700",
  derived: "text-fg-muted",
  predicted: "text-fg-muted",
  synthetic: "text-red-700",
};

/**
 * Everything a reader needs to decide whether to believe one number.
 *
 * The reviewer's list, and each line is here because its absence was doing
 * work: without units a Tm and a pH_T are both "a number"; without the domain
 * an extrapolation reads like an answer; without the training-row count a head
 * fitted to 160 measurements in one buffer looks like one fitted to 2,274
 * across 261.
 */
function HeadCard({ name, head }: { name: string; head: HeadInfo }) {
  const inert = head.inert_axes ?? [];
  const responds = head.responds_to ?? [];
  return (
    <details
      className="rounded-[var(--radius-md)] border border-border bg-bg/40 p-2"
      data-testid="head-card"
      data-head={name}
    >
      <summary className="flex cursor-pointer items-center justify-between gap-2 text-xs">
        <span className="font-mono text-fg">{name}</span>
        <span className={`text-[0.65rem] ${SEMANTICS_TONE[head.target_semantics] ?? ""}`}>
          {head.target_semantics_label ?? head.target_semantics}
        </span>
      </summary>

      <dl className="mt-2 grid grid-cols-[auto_1fr] gap-x-3 gap-y-1 text-[0.7rem]">
        <dt className="text-fg-subtle">Predicts</dt>
        <dd className="text-fg">
          {head.target}
          {head.units ? ` (${head.units})` : ""} · {head.task}
        </dd>

        <dt className="text-fg-subtle">Training rows</dt>
        <dd className="text-fg">{head.n_training_rows ?? "—"}</dd>

        {/* Two different statements, kept apart. */}
        <dt className="text-fg-subtle">Responds to</dt>
        <dd className="text-fg">
          {responds.length ? responds.join(", ") : "no condition axis"}
        </dd>

        {inert.length ? (
          <>
            <dt className="text-fg-subtle">No response to</dt>
            <dd className="text-fg-muted">{inert.join(", ")}</dd>
          </>
        ) : null}

        {head.predicted_axes?.length ? (
          <>
            <dt className="text-fg-subtle">Outputs</dt>
            <dd className="text-fg-muted">
              {head.predicted_axes.join(", ")} — an output, not an input
            </dd>
          </>
        ) : null}

        <dt className="text-fg-subtle">Calibration</dt>
        <dd className="text-fg-muted">{head.calibration_scope ?? "—"}</dd>
      </dl>

      {head.claim ? (
        <p className="mt-2 text-[0.68rem] leading-relaxed text-fg-subtle">{head.claim}</p>
      ) : null}
    </details>
  );
}

export function Workbench() {
  const backend = useNA((s) => s.backend);
  const [meta, setMeta] = useState<ServiceInfo | null>(null);
  const [tab, setTab] = useState<string>("sequence_evidence");

  useEffect(() => {
    if (backend !== "connected") return;
    let live = true;
    info()
      .then((m) => {
        if (live) setMeta(m as ServiceInfo);
      })
      .catch(() => {
        if (live) setMeta(null);
      });
    return () => {
      live = false;
    };
  }, [backend]);

  if (backend !== "connected" || !meta) return null;

  const current = meta.workspaces?.find((w) => w.id === tab) ?? meta.workspaces?.[0];

  return (
    <div className="space-y-3" data-testid="workbench">
      {/* Rides above everything, because a page of numbers from a model trained
          on generated data must not be readable without seeing this. */}
      {meta.synthetic_model ? (
        <div
          className="rounded-[var(--radius-md)] border border-red-500/50 bg-red-500/10 p-3"
          data-testid="synthetic-banner"
        >
          <p className="text-xs font-semibold text-red-700">
            {meta.synthetic_model.warning}
          </p>
          <p className="mt-1 text-[0.7rem] text-red-700/80">{meta.synthetic_model.note}</p>
        </div>
      ) : null}

      <div className="px-1 pt-1 pb-2">
        <p className="text-xs font-medium text-accent">EXPLORATION TOOLS</p>
        <h2 className="mt-1 text-xl font-semibold text-fg">What would you like to compare?</h2>
        <p className="mt-1 text-sm text-fg-muted">Choose a workflow. Each uses the sequence and buffer you set.</p>
      </div>
      <Tabs defaultValue="variants" className="min-w-0">
        <TabsList aria-label="Comparison workflows" className="mb-4 flex h-auto flex-wrap justify-start gap-1 bg-surface-2 p-1.5">
          <TabsTrigger value="variants">Variant join</TabsTrigger>
          <TabsTrigger value="mutations">Mutations</TabsTrigger>
          <TabsTrigger value="conditions">Conditions</TabsTrigger>
          <TabsTrigger value="batch">Batch</TabsTrigger>
          <TabsTrigger value="models">Model guide</TabsTrigger>
        </TabsList>
        <TabsContent forceMount value="variants" className="data-[state=inactive]:hidden"><VariantPanel /></TabsContent>
        <TabsContent forceMount value="mutations" className="data-[state=inactive]:hidden"><MutationExplorer /></TabsContent>
        <TabsContent forceMount value="conditions" className="data-[state=inactive]:hidden"><ConditionResponse /></TabsContent>
        <TabsContent forceMount value="batch" className="data-[state=inactive]:hidden"><BatchPanel /></TabsContent>
        <TabsContent forceMount value="models" className="data-[state=inactive]:hidden">

      <section className="rounded-[var(--radius-lg)] border border-border bg-surface p-4">
        <header>
          <h3 className="flex items-center gap-2 text-sm font-medium text-fg">
            <Info className="h-4 w-4" aria-hidden />
            Model guide
          </h3>
          <p className="text-xs text-fg-subtle">{meta.headline}</p>
        </header>

        <div
          className="mt-3 flex flex-wrap gap-1.5"
          role="tablist"
          aria-label="Head groups"
        >
          {meta.workspaces?.map((w) => (
            <button
              key={w.id}
              type="button"
              role="tab"
              aria-selected={tab === w.id}
              onClick={() => setTab(w.id)}
              data-testid="workspace-tab"
              data-workspace={w.id}
              className={`rounded-full border px-2.5 py-1 text-xs transition ${
                tab === w.id
                  ? "border-accent bg-accent/15 text-fg"
                  : "border-border text-fg-muted hover:text-fg"
              }`}
            >
              {w.label}
              <span className="ml-1 text-fg-subtle">{w.heads.length}</span>
            </button>
          ))}
        </div>

        {current ? (
          <div className="mt-3">
            <p className="text-[0.7rem] leading-relaxed text-fg-subtle">{current.note}</p>
            <div className="mt-2 space-y-1.5">
              {current.heads.length ? (
                current.heads.map((name) =>
                  meta.heads[name] ? (
                    <HeadCard key={name} name={name} head={meta.heads[name]} />
                  ) : null,
                )
              ) : (
                <p className="text-xs text-fg-muted">
                  No head in this model falls into this group.
                </p>
              )}
            </div>
          </div>
        ) : null}

        <p className="mt-3 flex gap-2 text-[0.65rem] leading-relaxed text-fg-subtle">
          <ChevronRight className="mt-0.5 h-3 w-3 shrink-0" aria-hidden />
          <span>
            Groups and per-head capabilities are read from the service, not
            listed here. Model {meta.model_version}, schema{" "}
            {meta.prediction_schema_version ?? "—"}.
          </span>
        </p>
      </section>
        </TabsContent>
      </Tabs>
    </div>
  );
}
