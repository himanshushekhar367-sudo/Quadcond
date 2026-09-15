import { HelpDialog } from "@/components/panels/HelpDialog";
import { useEffect } from "react";
import { LeftPanel, RightPanel } from "@/components/panels/Controls";
import { ViewerPane } from "@/components/viewer/ViewerPane";
import { Workbench } from "@/components/panels/Workbench";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { useNA } from "@/lib/na/store";

export function AppShell() {
  const polymer = useNA((s) => s.polymer);
  const n = useNA((s) => s.result.sequence.length);
  const backend = useNA((s) => s.backend);
  const bootstrap = useNA((s) => s.bootstrap);

  useEffect(() => {
    bootstrap();
  }, [bootstrap]);

  return (
    <div className="flex h-dvh overflow-hidden flex-col text-fg relative bg-transparent">
      
      <header className="glass relative z-10 gap-3 px-4 py-3 mx-4 mt-4 rounded-2xl flex items-center justify-between border border-border shadow-lg">
        <div className="flex flex-col">
          <p className="text-[0.65rem] tracking-[0.2em] text-accent font-semibold uppercase">QuadCond workbench</p>
          <h1 className="font-display text-2xl font-bold tracking-tight text-fg drop-shadow-md">AENNA<span className="text-accent">3D</span></h1>
        </div>
        
        <div className="hidden md:block flex-1 max-w-xl mx-8">
          <p className="text-xs leading-relaxed text-fg-muted">
            A viewer for non-canonical DNA structures. One sequence, its candidate
            folds drawn at published helical parameters, each annotated with what
            QuadCond predicts for it under the buffer you set.
          </p>
        </div>

        <div className="hidden sm:flex items-center gap-2 bg-surface-2 px-4 py-2 rounded-lg border border-border">
          <div className={`h-2 w-2 rounded-full ${backend === "connected" ? "bg-emerald-400" : "bg-amber-400"}`} title={backend === "connected" ? "Prediction service connected" : "Prediction service unavailable"} />
          <p className="font-mono text-xs tabular-nums text-fg-subtle font-medium">
            {polymer} · <span className="text-fg">{n} nt</span>
          </p>
        </div>
        <HelpDialog />
      </header>

      {/* Desktop Layout */}
      <div className="hidden min-h-0 flex-1 lg:grid lg:grid-cols-[320px_1fr_340px] gap-4 p-4 relative z-10">
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
            <div className="border-b border-border p-2">
              <TabsList className="w-full bg-surface-3">
                <TabsTrigger
                  className="flex-1 data-[state=active]:bg-accent data-[state=active]:text-white"
                  value="structure"
                >
                  Structure
                </TabsTrigger>
                <TabsTrigger
                  className="flex-1 data-[state=active]:bg-accent data-[state=active]:text-white"
                  value="workbench"
                  data-testid="workbench-tab"
                >
                  Compare
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
              {/* Subtle inner glow for the 3D viewer container */}
              <div className="absolute inset-0 border border-white/5 rounded-2xl pointer-events-none group-hover:border-accent/20 transition-colors duration-500" />
            </TabsContent>
            <TabsContent value="workbench" className="min-h-0 flex-1 overflow-y-auto p-3">
              <Workbench />
            </TabsContent>
          </Tabs>
        </main>
        
        <aside className="glass-panel flex flex-col overflow-hidden">
          <RightPanel />
        </aside>
      </div>

      {/* Mobile Layout */}
      <div className="flex min-h-0 flex-1 flex-col lg:hidden p-4 relative z-10">
        <Tabs defaultValue="view" className="glass-panel flex min-h-0 flex-1 flex-col overflow-hidden">
          <div className="border-b border-border p-2">
            <TabsList className="w-full bg-surface-3">
              <TabsTrigger className="flex-1 data-[state=active]:bg-accent data-[state=active]:text-white" value="seq">
                Sequence
              </TabsTrigger>
              <TabsTrigger className="flex-1 data-[state=active]:bg-accent data-[state=active]:text-white" value="view">
                3D
              </TabsTrigger>
              <TabsTrigger className="flex-1 data-[state=active]:bg-accent data-[state=active]:text-white" value="ens">
                Evidence
              </TabsTrigger>
              <TabsTrigger
                className="flex-1 data-[state=active]:bg-accent data-[state=active]:text-white"
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
