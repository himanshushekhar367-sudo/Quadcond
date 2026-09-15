import { useEffect, useState, type ComponentType } from "react";
import { useNA } from "@/lib/na/store";
import { KIND_LABEL } from "@/lib/na/types";
import type { GeometryOnlyMember } from "@/lib/na/geometry-only";

type SceneProps = {
  members: GeometryOnlyMember[];
  overlay: boolean;
  autoRotate: boolean;
};

export function ViewerPane() {
  const [Scene, setScene] = useState<ComponentType<SceneProps> | null>(null);
  const result = useNA((s) => s.result);
  const selectedId = useNA((s) => s.selectedId);
  const overlay = useNA((s) => s.overlay);
  const autoRotate = useNA((s) => s.autoRotate);

  useEffect(() => {
    let live = true;
    import("./NaScene").then((m) => {
      if (live) setScene(() => m.NaScene);
    });
    return () => {
      live = false;
    };
  }, []);

  // `result.members`, not `result.ensemble`. The store stopped producing a
  // local ensemble when the viewer's numeric layer was removed, this component
  // was never updated, and `[]` renders as "No structure" -- so every fold
  // vanished from the 3D pane, connected and offline alike, while the panels
  // beside it showed correct predictions. A shape change on one side of a
  // boundary with no type error on the other.
  // Exact id, no fallback. `?? members[0]` is what turned every unmatched
  // selection into "whatever happens to be first" — which for a G-rich sequence
  // is the parallel G4, so clicking the i-motif card drew a quadruplex. An
  // unmatched id now shows nothing, and the panel says why.
  const selected = result.members.find((m) => m.id === selectedId);
  const members = overlay ? result.members : selected ? [selected] : [];

  return (
    <div className="relative h-full min-h-[320px] overflow-hidden rounded-[var(--radius-xl)] border border-border bg-bg">
      {Scene ? (
        <Scene members={members} overlay={overlay} autoRotate={autoRotate} />
      ) : (
        <div className="flex h-full items-center justify-center text-sm text-fg-muted">
          Loading viewer
        </div>
      )}
      <div className="pointer-events-none absolute inset-x-0 top-0 flex items-start justify-between p-4">
        <div className="rounded-[var(--radius-md)] border border-border bg-bg/80 px-3 py-2">
          <p className="text-[0.65rem] tracking-wide text-fg-subtle uppercase">Active fold</p>
          {/* Kind alone does not identify a member. Two members of the same
              kind differ by strand and topology, and the id-namespace defect
              this attribute exists for produced exactly that pair: the right
              kind on the wrong strand. The gate asserts all three. */}
          <p
            className="text-sm text-fg"
            data-testid="active-fold"
            data-kind={selected?.kind ?? ""}
            data-strand={selected?.strand ?? ""}
            data-topology={selected?.topology ?? ""}
            data-member-id={selected?.id ?? ""}
          >
            {selected ? KIND_LABEL[selected.kind] : "—"}
          </p>
          {selected ? (
            <p className="text-[0.6rem] text-fg-subtle">
              schematic archetype · strand {selected.strand}
            </p>
          ) : null}
        </div>
        <div className="rounded-[var(--radius-md)] border border-border bg-bg/80 px-3 py-2 text-right">
          <p className="text-[0.65rem] tracking-wide text-fg-subtle uppercase">Bases</p>
          <p className="font-mono text-xs text-fg">
            <span className="text-base-a">A</span>{" "}
            <span className="text-base-t">T/U</span>{" "}
            <span className="text-base-g">G</span>{" "}
            <span className="text-base-c">C</span>
          </p>
        </div>
      </div>
      <p className="pointer-events-none absolute bottom-3 left-4 text-[0.65rem] text-fg-subtle">
        Drag to orbit · scroll to zoom
      </p>
    </div>
  );
}
