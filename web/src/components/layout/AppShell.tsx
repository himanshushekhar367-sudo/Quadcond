import { HelpDialog } from "@/components/panels/HelpDialog";
import { useEffect } from "react";
import { LeftPanel, RightPanel } from "@/components/panels/Controls";
import { ViewerPane } from "@/components/viewer/ViewerPane";
import { Workbench } from "@/components/panels/Workbench";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { useNA } from "@/lib/na/store";
import { Dna, FlaskConical, ArrowRight } from "lucide-react";

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
      
      <header className="app-header relative z-10 flex items-center justify-between gap-3 border-b border-border bg-surface px-5 py-4">
        <div className="flex items-center gap-3">
          <div className="flex h-11 w-11 items-center justify-center rounded-xl bg-accent/10 text-accent"><Dna size={26} aria-hidden /></div>
          <div>
            <h1 className="font-display text-xl font-semibold tracking-tight text-fg">QuadCond <span className="font-normal text-fg-subtle">/</span> <span className="text-accent">AENNA-3D</span></h1>
            <p className="mt-0.5 text-xs text-fg-muted">DNA structure &amp; variant evidence</p>
          </div>
        </div>
        
        <div className="hidden xl:flex items-center gap-3 text-xs text-fg-muted" aria-label="Workflow">
          <span>1. Add a sequence</span><ArrowRight size={13} aria-hidden />
          <span>2. Set conditions</span><ArrowRight size={13} aria-hidden />
          <span>3. Explore evidence</span>
        </div>

        <div className="hidden sm:flex items-center gap-2 bg-surface-2 px-3 py-2 rounded-full border border-border">
          <div className={`h-2 w-2 rounded-full ${backend === "connected" ? "bg-emerald-400" : "bg-amber-400"}`} title={backend === "connected" ? "Prediction service connected" : "Prediction service unavailable"} />
          <p className="text-xs tabular-nums text-fg-muted">
            {backend === "connected" ? "Model connected" : "Model unavailable"} <span className="ml-2 font-mono">{polymer} · {n} nt</span>
          </p>
        </div>
        <HelpDialog />
      </header>

      {/* Desktop Layout */}
      <div className="hidden min-h-0 flex-1 lg:grid lg:grid-cols-[280px_minmax(0,1fr)_320px] xl:grid-cols-[296px_minmax(0,1fr)_340px] gap-4 p-4 relative z-10">
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
            <div className="border-b border-border p-3">
              <TabsList className="w-full bg-surface-3">
                <TabsTrigger
                  className="flex-1 data-[state=active]:bg-accent data-[state=active]:text-white"
                  value="structure"
                >
                  <Dna size={15} className="mr-2" aria-hidden /> Structure
                </TabsTrigger>
                <TabsTrigger
                  className="flex-1 data-[state=active]:bg-accent data-[state=active]:text-white"
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
              {/* Subtle inner glow for the 3D viewer container */}
              <div className="absolute inset-0 border border-border rounded-2xl pointer-events-none group-hover:border-accent/20 transition-colors duration-500" />
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
