import { EXAMPLES } from "@/lib/na/examples";
import { DEFAULT_CONDITIONS, KIND_LABEL, type Conditions } from "@/lib/na/types";
import { useNA } from "@/lib/na/store";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { ScrollArea } from "@/components/ui/scroll-area";
import { Slider } from "@/components/ui/slider";
import { Separator } from "@/components/ui/separator";
import { 
  RotateCcw, Thermometer, Zap, Droplets, FlaskConical,
  Dna, Activity, Library, Layers, Sparkles, ChevronRight, RefreshCw
} from "lucide-react";
import {
  GenomicEvidencePanel,
  LocusOverlapPanel,
  EvidenceCards,
  RefusalNotice,
  StructureProvenance,
  WithheldFolds,
} from "./EvidencePanels";

function ConditionRow({
  label,
  unit,
  value,
  min,
  max,
  step,
  onChange,
  icon: Icon,
  hint,
}: {
  label: string;
  unit: string;
  value: number;
  min: number;
  max: number;
  step: number;
  onChange: (v: number) => void;
  icon: any;
  hint?: string;
}) {
  return (
    <div className="space-y-3 p-3 rounded-xl bg-black/20 border border-white/5 hover:border-accent/30 hover:bg-black/40 transition-all duration-300 group">
      <div className="flex items-center justify-between gap-3">
        <div className="flex items-center gap-2 text-fg-muted group-hover:text-accent transition-colors">
          <Icon className="w-4 h-4 drop-shadow-md" />
          <span className="text-xs font-semibold uppercase tracking-wider">{label}</span>
        </div>
        <div className="px-2 py-0.5 rounded-md bg-black/40 border border-white/10 shadow-inner">
          <span className="font-mono text-xs tabular-nums text-accent font-bold">
            {value}
            <span className="text-fg-subtle text-[10px] ml-1">{unit}</span>
          </span>
        </div>
      </div>
      <Slider
        min={min}
        max={max}
        step={step}
        value={[value]}
        onValueChange={(v) => onChange(v[0] ?? value)}
        className="cursor-pointer"
      />
      {hint ? (
        <p className="text-[10px] leading-snug text-fg-subtle">{hint}</p>
      ) : null}
    </div>
  );
}

export function LeftPanel() {
  const sequence = useNA((s) => s.sequence);
  const polymer = useNA((s) => s.polymer);
  const conditions = useNA((s) => s.conditions);
  const setSequence = useNA((s) => s.setSequence);
  const setPolymer = useNA((s) => s.setPolymer);
  const setCondition = useNA((s) => s.setCondition);
  const loadExample = useNA((s) => s.loadExample);

  return (
    <ScrollArea className="h-full">
      <div className="space-y-8 p-5">
        
        {/* Sequence Input Section */}
        <div className="space-y-3">
          <div className="flex items-center gap-2 text-fg-subtle">
            <Dna className="w-4 h-4 text-accent" />
            <p className="text-xs font-bold tracking-widest uppercase">Sequence Input</p>
          </div>
          
          <div className="flex gap-2 p-1 rounded-lg bg-black/30 border border-white/5">
            {(["DNA", "RNA"] as const).map((p) => (
              <Button
                key={p}
                size="sm"
                variant="ghost"
                className={`flex-1 text-xs font-bold tracking-wider transition-all ${
                  polymer === p 
                    ? "bg-accent/20 text-accent border border-accent/50 shadow-[0_0_15px_rgba(59,130,246,0.3)]" 
                    : "text-fg-muted hover:text-fg hover:bg-white/5"
                }`}
                onClick={() => setPolymer(p)}
              >
                {p}
              </Button>
            ))}
          </div>
          
          <div className="relative group">
            <div className="absolute -inset-0.5 bg-gradient-to-r from-accent/0 to-accent/0 group-focus-within:from-accent/40 group-focus-within:to-purple-500/40 rounded-xl blur opacity-50 transition duration-500"></div>
            <textarea
              value={sequence}
              onChange={(e) => setSequence(e.target.value)}
              spellCheck={false}
              rows={4}
              className="relative w-full resize-y rounded-xl border border-white/10 bg-black/60 px-4 py-3 font-mono text-sm leading-relaxed text-fg outline-none focus:border-accent/50 backdrop-blur-md transition-all shadow-inner"
              aria-label="Nucleic acid sequence"
            />
          </div>
          <div className="flex justify-end">
            <Badge variant="outline" className="bg-black/50 border-white/10 font-mono text-[10px] text-accent">
              {sequence.replace(/[^ACGTUacgtu]/g, "").length} nt
            </Badge>
          </div>
        </div>

        <Separator className="bg-white/10" />

        {/* Library Section */}
        <details className="group [&_summary::-webkit-details-marker]:hidden">
          <summary className="flex items-center justify-between cursor-pointer list-none text-fg-subtle hover:text-fg transition-colors select-none">
            <div className="flex items-center gap-2">
              <Library className="w-4 h-4 text-purple-400" />
              <p className="text-xs font-bold tracking-widest uppercase">Presets Library</p>
            </div>
            <ChevronRight className="w-4 h-4 transition-transform group-open:rotate-90" />
          </summary>
          
          <div className="grid grid-cols-1 gap-3 mt-4">
            {EXAMPLES.map((ex) => (
              <button
                key={ex.id}
                type="button"
                onClick={() => loadExample(ex.id)}
                className="group/btn relative overflow-hidden rounded-xl border border-white/10 bg-gradient-to-b from-surface-2 to-black/60 p-3 text-left transition-all duration-300 hover:border-purple-400/50 hover:shadow-[0_0_20px_rgba(168,85,247,0.15)]"
              >
                <div className="absolute inset-0 bg-gradient-to-r from-purple-500/0 via-purple-500/0 to-purple-500/0 group-hover/btn:via-purple-500/10 transition-all duration-700 translate-x-[-100%] group-hover/btn:translate-x-[100%]"></div>
                <div className="flex items-center justify-between gap-2 relative z-10">
                  <span className="text-sm font-semibold text-fg group-hover/btn:text-purple-300 transition-colors">{ex.name}</span>
                  <div className="flex items-center gap-1">
                    {ex.coverage && ex.coverage !== "full" ? (
                      <Badge
                        variant="outline"
                        className={`text-[9px] ${
                          ex.coverage === "refused"
                            ? "border-red-500/40 bg-red-500/10 text-red-300"
                            : "border-amber-500/40 bg-amber-500/10 text-amber-300"
                        }`}
                      >
                        {ex.coverage === "refused"
                          ? "refused"
                          : ex.coverage === "partial"
                            ? "partly covered"
                            : "not modelled"}
                      </Badge>
                    ) : null}
                    <Badge className="bg-black/50 border-purple-500/30 text-purple-200 text-[10px]">{ex.polymer}</Badge>
                  </div>
                </div>
                <p className="mt-1.5 text-xs text-fg-muted relative z-10">{ex.note}</p>
              </button>
            ))}
          </div>
        </details>

        <Separator className="bg-white/10" />

        {/* Conditions Section */}
        <details className="group [&_summary::-webkit-details-marker]:hidden">
          <summary className="flex items-center justify-between cursor-pointer list-none text-fg-subtle hover:text-fg transition-colors select-none">
            <div className="flex items-center gap-2">
              <Activity className="w-4 h-4 text-green-400" />
              <p className="text-xs font-bold tracking-widest uppercase">
                Buffer conditions
              </p>
            </div>
            <ChevronRight className="w-4 h-4 transition-transform group-open:rotate-90" />
          </summary>
          
          <div className="space-y-2 mt-4">
            <ConditionRow
              icon={Droplets}
              label="pH"
              unit=""
              value={round(conditions.ph, 2)}
              min={4}
              max={8.5}
              step={0.1}
              onChange={(v) => setCondition("ph", v)}
            />

            {/*
              K+ and Na+ are separate inputs because they are separate
              chemistry. On the same sequence at the same concentration,
              potassium melts a quadruplex about 13 degC higher than sodium --
              measured across 411 sequence-matched pairs, and reproduced
              out-of-fold by the g4_tm head at 94% sign agreement. A single
              "monovalent" slider mapped to K+ answered a potassium question
              for a user who typed a sodium buffer.
            */}
            <ConditionRow
              icon={Zap}
              label="Potassium K+"
              unit="mM"
              value={conditions.k}
              min={0}
              max={500}
              step={1}
              onChange={(v) => setCondition("k", v)}
            />
            <ConditionRow
              icon={Zap}
              label="Sodium Na+"
              unit="mM"
              value={conditions.na}
              min={0}
              max={500}
              step={1}
              onChange={(v) => setCondition("na", v)}
            />
            <ConditionRow
              icon={Zap}
              label="Li+ / NH4+"
              unit="mM"
              value={conditions.li_nh4}
              min={0}
              max={500}
              step={1}
              onChange={(v) => setCondition("li_nh4", v)}
              hint="Lithium is the standard destabilising control: it does not fit the tetrad channel."
            />
            <ConditionRow
              icon={Sparkles}
              label="Magnesium Mg2+"
              unit="mM"
              value={round(conditions.mg, 2)}
              min={0}
              max={20}
              step={0.1}
              onChange={(v) => setCondition("mg", v)}
            />
            <ConditionRow
              icon={FlaskConical}
              label="Crowder"
              unit="% w/v"
              value={conditions.crowder_pct}
              min={0}
              max={40}
              step={1}
              onChange={(v) => setCondition("crowder_pct", v)}
            />
            <ConditionRow
              icon={Layers}
              label="Strand conc."
              unit="µM"
              value={round(conditions.strand_conc, 1)}
              min={0.1}
              max={100}
              step={0.1}
              onChange={(v) => setCondition("strand_conc", v)}
              hint="Matters for anything intermolecular; the atlas records it because published melting curves depend on it."
            />
            <ConditionRow
              icon={Thermometer}
              label="Temperature"
              unit="°C"
              value={conditions.temperature}
              min={4}
              max={95}
              step={1}
              onChange={(v) => setCondition("temperature", v)}
            />
          </div>

          <Button
            variant="outline"
            className="w-full mt-4 bg-transparent border-white/10 hover:bg-white/5 text-xs tracking-wider uppercase text-fg-muted hover:text-fg transition-all"
            onClick={() => {
              (Object.keys(DEFAULT_CONDITIONS) as (keyof Conditions)[]).forEach((k) =>
                setCondition(k, DEFAULT_CONDITIONS[k]),
              );
            }}
          >
            <RotateCcw className="w-3 h-3 mr-2" />
            Reset to nominal cytosol
          </Button>
          <p className="mt-2 text-[10px] leading-relaxed text-fg-subtle">
            Nominal, not measured. Every field here is a choice; QuadCond marks
            the ones you never set as imputed on the response.
          </p>
        </details>
      </div>
    </ScrollArea>
  );
}

export function RightPanel() {
  const result = useNA((s) => s.result);
  const selectedId = useNA((s) => s.selectedId);
  const selectStructure = useNA((s) => s.selectStructure);
  const overlay = useNA((s) => s.overlay);
  const setOverlay = useNA((s) => s.setOverlay);
  const autoRotate = useNA((s) => s.autoRotate);
  const setAutoRotate = useNA((s) => s.setAutoRotate);
  const visMode = useNA((s) => s.visMode);
  const setVisMode = useNA((s) => s.setVisMode);
  const backend = useNA((s) => s.backend);
  const backendMessage = useNA((s) => s.backendMessage);
  const busy = useNA((s) => s.busy);

  const quadcond = useNA((s) => s.quadcond);
  const qcCards = quadcond?.evidence?.length ? quadcond.evidence : null;
  const joint = quadcond?.locus_peak_overlap_state ?? null;
  // `result.members`. This read was `result.ensemble` — always `[]` once the
  // local ensemble was removed — so `selected` was permanently undefined and
  // the structure-provenance card and dot-bracket notation never rendered at
  // all. Deleting `ensemble` from GeometryResult is what surfaced it.
  const selected = result.members.find((m) => m.id === selectedId);

  return (
    <ScrollArea className="h-full">
      <div className="space-y-8 p-5">
        
        {/*
          What the backend is actually doing. The panel this replaces showed
          "QuaDB ML Model Inferring..." on a 600 ms setTimeout and "ML Model
          Ready" the rest of the time -- there was no model and no inference.
          There is now, so the indicator reports its real state, including the
          state where it is unreachable and the numbers on screen are local.
        */}
        <div
          className={`flex items-center justify-between gap-2 p-2 rounded-lg border transition-all duration-300 ${
            backend === "connected"
              ? "bg-black/30 border-white/5"
              : backend === "offline"
                ? "bg-amber-500/10 border-amber-500/40"
                : "bg-black/30 border-white/5"
          }`}
        >
          <div className="flex items-center gap-2">
            <Activity
              className={`w-4 h-4 ${busy ? "text-accent animate-pulse" : backend === "offline" ? "text-amber-400" : "text-fg-subtle"}`}
            />
            <span className="text-[10px] font-bold tracking-widest uppercase text-fg-subtle">
              {busy
                ? "predicting…"
                : backend === "connected"
                  ? "QuadCond connected"
                  : backend === "offline"
                    ? "QuadCond unreachable — geometry only"
                    : "QuadCond not contacted"}
            </span>
          </div>
        </div>
        {backend === "offline" ? (
          <p className="-mt-6 text-[10px] leading-relaxed text-amber-300/80">
            The structures below are drawn at published helical parameters, but
            no calibrated number is available. Nothing is substituted for them.
            {backendMessage ? ` (${backendMessage})` : null}
          </p>
        ) : null}

        {/* Viewer Controls */}
        <div className="flex flex-col gap-2 p-2 rounded-lg bg-black/30 border border-white/5">
          <div className="flex items-center justify-between">
            <Button 
              size="sm" 
              variant="ghost" 
              className={`flex-1 text-xs font-semibold tracking-wider transition-all ${overlay ? "bg-accent/20 text-accent" : "text-fg-muted"}`}
              onClick={() => setOverlay(!overlay)}
            >
              <Layers className="w-3 h-3 mr-2" />
              Overlay Folds
            </Button>
            <div className="w-px h-4 bg-white/10 mx-1" />
            <Button
              size="sm"
              variant="ghost"
              className={`flex-1 text-xs font-semibold tracking-wider transition-all ${autoRotate ? "bg-accent/20 text-accent" : "text-fg-muted"}`}
              onClick={() => setAutoRotate(!autoRotate)}
            >
              <RefreshCw className={`w-3 h-3 mr-2 ${autoRotate ? 'animate-spin-slow' : ''}`} />
              Auto-Rotate
            </Button>
          </div>
          <div className="flex items-center gap-1 border-t border-white/5 pt-2">
            {(['ball-and-stick', 'space-filling', 'licorice'] as const).map(mode => (
              <Button 
                key={mode} 
                size="sm" 
                variant="ghost" 
                className={`flex-1 text-[9px] font-bold tracking-widest uppercase transition-all ${visMode === mode ? 'bg-purple-500/20 text-purple-300 border border-purple-500/30' : 'text-fg-muted hover:bg-white/5'}`}
                onClick={() => setVisMode(mode)}
              >
                {mode.replace(/-/g, ' ')}
              </Button>
            ))}
          </div>
        </div>

        {/* Structural evidence — not an ensemble; see EvidencePanels.tsx */}
        <details className="group [&_summary::-webkit-details-marker]:hidden" open>
          <summary className="flex cursor-pointer select-none items-center justify-between text-fg-subtle transition-colors hover:text-fg">
            <div className="flex items-center gap-2">
              <Activity className="w-4 h-4 text-accent" />
              <p className="text-xs font-bold tracking-widest uppercase">Structural evidence</p>
            </div>
            <ChevronRight className="w-4 h-4 transition-transform group-open:rotate-90" />
          </summary>

          {result.refusal ? (
            <div className="mt-4">
              <RefusalNotice reason={result.refusal} />
            </div>
          ) : qcCards ? (
            <>
              <EvidenceCards
                cards={qcCards}
                note={quadcond?.evidence_note}
                syntheticModel={quadcond?.synthetic_model}
              />
              {quadcond?.condition_imputed_fields?.length ? (
                <p className="mt-3 text-[10px] leading-relaxed text-fg-subtle">
                  Imputed (you did not set these):{" "}
                  {quadcond.condition_imputed_fields.join(", ")}. A number computed
                  under a defaulted buffer is not a prediction under a measured one.
                </p>
              ) : null}
            </>
          ) : (
            /*
              Offline: structures, and no numbers at all.

              What used to be here was a local ensemble with ΔG values and
              probability bars, under a banner reading "geometry only — nothing
              is substituted". The numbers were substituted, and they were bad
              in three separate ways: the energy terms read total monovalent
              concentration so 150 mM K+ and 150 mM Na+ gave a byte-identical
              answer; the ΔG calibration was never applied on this path, leaving
              a 46 kcal/mol span that saturated the Boltzmann clamp at
              0.5/0.5/0.000; and the coil reference came out at +3 kcal/mol
              rather than 0, so every value floated on an unstated offset.

              A viewer with no backend has structures to draw and nothing to
              say about stability. That is the honest state, and it is now the
              rendered one.
            */
            <div className="mt-4 space-y-3">
              <div className="rounded-lg border border-amber-500/30 bg-amber-500/5 p-3">
                <p className="text-[11px] leading-relaxed text-amber-200/90">
                  {backend === "connected" && !busy ? (
                    <>QuadCond is connected, but returned no structural evidence
                    cards for this sequence under its implemented rules. This is
                    an applicability result, not evidence that the DNA cannot form
                    a structure. The drawings below are schematic.</>
                  ) : busy ? (
                    <>Waiting for QuadCond predictions. The drawings below are schematic.</>
                  ) : (<>
                  QuadCond is not reachable, so there are no predicted values —
                  no melting temperature, no transitional pH, no probabilities.
                  The structures below are drawn at published helical parameters
                  from the motifs in your sequence. Nothing has been substituted
                  for the missing numbers.
                  </>)}
                </p>
              </div>
              <ul className="space-y-2">
                {result.members.map((m) => {
                  const active = m.id === selectedId;
                  return (
                    <li key={m.id}>
                      <button
                        type="button"
                        onClick={() => selectStructure(m.id)}
                        className={`w-full rounded-xl border p-3 text-left transition-all ${
                          active
                            ? "border-accent/50 bg-accent/5"
                            : "border-white/5 bg-black/40 hover:border-white/20 hover:bg-black/60"
                        }`}
                      >
                        <div className="flex items-center justify-between gap-2">
                          <span className={`text-sm font-bold ${active ? "text-accent" : "text-fg"}`}>
                            {KIND_LABEL[m.kind]}
                          </span>
                          <div className="flex items-center gap-1">
                            <span className="font-mono text-[9px] text-fg-subtle">
                              strand {m.strand}
                            </span>
                            <Badge
                              variant="outline"
                              className="border-white/10 text-[9px] tracking-wide text-fg-subtle"
                            >
                              {m.topology}
                            </Badge>
                          </div>
                        </div>
                        <p className="mt-2 text-[10px] leading-relaxed text-fg-subtle">
                          {m.description}
                        </p>
                      </button>
                    </li>
                  );
                })}
              </ul>
            </div>
          )}
        </details>

        {selected && (
          <>
            {/* Provenance is looked up for the strand actually drawn — a
                template match is a claim about a specific molecule, and the
                telomeric 22-mer and its complement are different ones. */}
            <StructureProvenance sequence={selected.strandSequence} kind={selected.kind} />
            <div className="rounded-xl border border-white/5 bg-black/30 p-4">
              <div className="mb-2 flex items-center justify-between gap-2">
                <p className="text-xs font-bold uppercase tracking-widest text-fg-subtle">
                  Dot-bracket notation
                </p>
                <span className="font-mono text-[10px] text-fg-subtle">
                  strand {selected.strand}
                </span>
              </div>
              <p className="break-all font-mono text-[11px] leading-relaxed text-accent">
                {selected.notation}
              </p>
              <p className="mt-2 break-all font-mono text-[10px] leading-relaxed text-fg-subtle">
                {selected.strandSequence}
              </p>
              <p className="mt-3 border-l-2 border-white/10 pl-2 text-[11px] text-fg-muted">
                {selected.description}
              </p>
            </div>
          </>
        )}

        <Separator className="bg-white/10" />

        {/*
          The protein-binder panel is gone rather than hidden. Its records were
          unreliable -- BLM (UniProt P54132) was stored as 249 residues against a
          true length of 1,417 -- and "% Match" was computed from motif class,
          not from any binding measurement. See src/lib/na/capabilities.ts.
        */}
        <GenomicEvidencePanel
          items={quadcond?.genomic_evidence ?? []}
          note={quadcond?.genomic_evidence_note}
        />

        {joint ? <LocusOverlapPanel joint={joint} /> : null}

        <Separator className="bg-white/10" />

        <WithheldFolds />

      </div>
    </ScrollArea>
  );
}

function round(n: number, d: number) {
  const p = 10 ** d;
  return Math.round(n * p) / p;
}
