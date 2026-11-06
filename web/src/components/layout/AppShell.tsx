import { HelpDialog } from "@/components/panels/HelpDialog";
import { useEffect } from "react";
import { LeftPanel, RightPanel } from "@/components/panels/Controls";
import { ViewerPane } from "@/components/viewer/ViewerPane";
import { Workbench } from "@/components/panels/Workbench";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { useNA } from "@/lib/na/store";
import { Dna, FlaskConical, ChevronRight, ShieldCheck } from "lucide-react";

export function AppShell() {
  const polymer = useNA((s) => s.polymer);
  const n = useNA((s) => s.result.sequence.length);
  const backend = useNA((s) => s.backend);
  const bootstrap = useNA((s) => s.bootstrap);
  const backendMessage = useNA((s) => s.backendMessage);
  const version = useNA((s) => s.quadcond?.quadcond_version);

  useEffect(() => {
    bootstrap();
  }, [bootstrap]);

  return (
    <div className="flex h-dvh overflow-hidden flex-col text-fg relative bg-transparent">
      
      <header className="app-header relative z-10 flex h-14 shrink-0 items-center justify-between gap-4 border-b border-border bg-surface px-4 shadow-[0_1px_2px_rgb(16_28_38/0.04)]">
        <div className="flex min-w-0 items-center gap-3">
          <div className="flex size-8 shrink-0 items-center justify-center rounded-md bg-accent text-white shadow-[0_1px_2px_rgb(16_28_38/0.2)]"><Dna size={18} aria-hidden /></div>
          <div className="min-w-0">
            <h1 className="font-display text-[15px] font-semibold leading-tight tracking-tight text-fg">
              <span className="font-medium text-fg-muted">QuadCond</span> <span className="font-normal text-border-strong">/</span> AENNA-3D
            </h1>
            <p className="truncate text-[11px] leading-tight text-fg-subtle">G-quadruplex &amp; i-motif evidence workbench</p>
          </div>
          {version ? <span className="hidden rounded-[var(--radius-xs)] border border-border px-1.5 py-0.5 font-mono text-[10px] text-fg-subtle md:inline">v{version}</span> : null}
        </div>

        <ol className="hidden items-center gap-2 text-xs text-fg-muted xl:flex" aria-label="Workflow">
          {["Add a sequence", "Set conditions", "Explore evidence"].map((step, i) => (
            <li key={step} className="flex items-center gap-2">
              {i > 0 ? <ChevronRight size={13} className="text-border-strong" aria-hidden /> : null}
              <span className="flex size-[18px] items-center justify-center rounded-full border border-border-strong bg-surface-2 font-mono text-[10px] font-medium text-fg-muted">{i + 1}</span>
              <span>{step}</span>
            </li>
          ))}
        </ol>

        <div className="flex items-center gap-2">
          <div className="hidden items-center overflow-hidden rounded-[var(--radius-md)] border border-border bg-surface-2 text-xs sm:flex">
            <span className="flex items-center gap-2 px-2.5 py-1.5">
              <span className={`size-2 rounded-full ${backend === "connected" ? "bg-ok shadow-[0_0_0_3px_rgb(31_122_77/0.15)]" : backend === "unknown" ? "bg-fg-subtle/60" : "bg-amber-500 shadow-[0_0_0_3px_rgb(245_158_11/0.18)]"}`} title={backend === "connected" ? "Prediction service connected" : backend === "unknown" ? "Contacting the prediction service" : "Prediction service unavailable"} />
              <span className="font-medium text-fg" data-testid="model-status">{backend === "connected" ? "Model connected" : backend === "unknown" ? "Connecting…" : "Model unavailable"}</span>
            </span>
            <span className="border-l border-border px-2.5 py-1.5 font-mono tabular-nums text-fg-muted">{polymer} · {n} nt</span>
          </div>
          <HelpDialog />
        </div>
      </header>

      {/* Desktop Layout */}
      <div className="hidden min-h-0 flex-1 lg:grid lg:grid-cols-[280px_minmax(0,1fr)_320px] xl:grid-cols-[296px_minmax(0,1fr)_340px] gap-3 p-3 relative z-10">
        <aside className="glass-panel flex flex-col overflow-hidden">
          <LeftPanel />
        </aside>
        
        {/*
          The centre is shared, not owned by the picture.

          The comparison tools first lived in the 340px right rail, where their
          tables clipped and their charts were unreadable — which is the same
          mistake in layout that the rest of this project keeps making in prose:
          the most persuasive element takes the space, and the evidence needed
          to act on it gets whatever is left. A schematic archetype is an
          explanation of a hypothesis; the mutation scan is the thing a
          researcher makes a decision from. They get equal room, and the tab
          strip means neither is buried under the other.
        */}
        <main className="glass-card flex min-w-0 flex-col overflow-hidden relative group">
          <Tabs defaultValue="structure" className="flex min-h-0 flex-1 flex-col">
            <div className="border-b border-border px-4">
              <TabsList className="h-11 w-full justify-start gap-6 rounded-none border-0 bg-transparent p-0">
                <TabsTrigger
                  className="h-11 rounded-none border-b-2 border-transparent px-0.5 text-[13px] data-[state=active]:border-accent data-[state=active]:bg-transparent data-[state=active]:text-fg data-[state=active]:shadow-none"
                  value="structure"
                >
                  <Dna size={15} className="mr-2" aria-hidden /> Structure
                </TabsTrigger>
                <TabsTrigger
                  className="h-11 rounded-none border-b-2 border-transparent px-0.5 text-[13px] data-[state=active]:border-accent data-[state=active]:bg-transparent data-[state=active]:text-fg data-[state=active]:shadow-none"
                  value="workbench"
                  data-testid="workbench-tab"
                >
                  <FlaskConical size={15} className="mr-2" aria-hidden /> Compare
                </TabsTrigger>
              </TabsList>
            </div>
            {/*
              `forceMount`, with the inactive state hidden by CSS rather than
              unmounted.

              Radix unmounts inactive content by default, which tore down the
              WebGL context on every switch to Compare and back. Two problems,
              one visible and one not: react-three-fiber threw
              "Cannot read properties of null (reading 'addEventListener')" as
              it detached from a node Radix had already removed, and rebuilding
              the context costs about four seconds on software rendering — so
              the tab that exists to stop the picture monopolising attention was
              punishing anyone who looked away from it.
            */}
            <TabsContent
              forceMount
              value="structure"
              className="relative min-h-0 flex-1 data-[state=inactive]:hidden"
            >
              <ViewerPane />
            </TabsContent>
            <TabsContent value="workbench" className="min-h-0 flex-1 overflow-y-auto bg-surface-2 p-3">
              <Workbench />
            </TabsContent>
          </Tabs>
        </main>
        
        <aside className="glass-panel flex flex-col overflow-hidden">
          <RightPanel />
        </aside>
      </div>

      {/* Status bar: the model identity is always on screen, not only in exports. */}
      <footer className="app-statusbar relative z-10 hidden h-7 shrink-0 items-center justify-between gap-4 border-t border-border bg-surface px-4 text-fg-subtle lg:flex">
        <span className="flex min-w-0 items-center gap-2 truncate">
          <span className={`size-1.5 shrink-0 rounded-full ${backend === "connected" ? "bg-ok" : backend === "unknown" ? "bg-fg-subtle/60" : "bg-amber-500"}`} aria-hidden />
          <span className="truncate font-mono">{backend === "connected" ? backendMessage : backend === "unknown" ? "Contacting the prediction service…" : "Prediction service not connected: schematic geometry only, no predicted values"}</span>
        </span>
        <span className="flex shrink-0 items-center gap-3">
          <span className="flex items-center gap-1.5"><ShieldCheck size={12} aria-hidden /> No accounts · requests are not stored</span>
          <span className="text-border-strong">|</span>
          <span>MIT licence</span>
        </span>
      </footer>

      {/* Mobile Layout */}
      <div className="flex min-h-0 flex-1 flex-col lg:hidden p-3 relative z-10">
        <Tabs defaultValue="view" className="glass-panel flex min-h-0 flex-1 flex-col overflow-hidden">
          <div className="border-b border-border p-2">
            <TabsList className="w-full">
              <TabsTrigger className="flex-1" value="seq">
                Sequence
              </TabsTrigger>
              <TabsTrigger className="flex-1" value="view">
                3D
              </TabsTrigger>
              <TabsTrigger className="flex-1" value="ens">
                Evidence
              </TabsTrigger>
              <TabsTrigger
                className="flex-1"
                value="compare"
              >
                Compare
              </TabsTrigger>
            </TabsList>
          </div>
          <TabsContent value="seq" className="min-h-0 flex-1 p-2">
            <div className="h-full">
              <LeftPanel />
            </div>
          </TabsContent>
          <TabsContent forceMount value="view" className="min-h-0 flex-1 relative data-[state=inactive]:hidden">
            <div className="h-full absolute inset-0">
              <ViewerPane />
            </div>
          </TabsContent>
          <TabsContent value="ens" className="min-h-0 flex-1 p-2">
            <div className="h-full">
              <RightPanel />
            </div>
          </TabsContent>
          <TabsContent value="compare" className="min-h-0 flex-1 overflow-y-auto p-2">
            <Workbench />
          </TabsContent>
        </Tabs>
      </div>
    </div>
  );
}
