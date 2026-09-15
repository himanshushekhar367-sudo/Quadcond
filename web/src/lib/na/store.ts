import { create } from "zustand";
import { EXAMPLES } from "./examples";
import { geometryOnly, memberForCard, type GeometryResult } from "./geometry-only";
import {
  DEFAULT_CONDITIONS,
  type Conditions,
  type Polymer,
} from "./types";
import {
  evidence as fetchEvidence,
  health as fetchHealth,
  type EvidenceResponse,
  type QuadcondConditions,
} from "../quadcond/client";

const first = EXAMPLES[0]!;

export type BackendStatus = "unknown" | "connected" | "offline";

export interface NAState {
  sequence: string;
  polymer: Polymer;
  conditions: Conditions;
  /**
   * Structures to draw. Local, because building coordinates from published
   * helical parameters is the viewer's own job — and it carries no numbers,
   * because computing those locally is not.
   */
  result: GeometryResult;
  /** Calibrated numbers from the QuadCond service. Null when it is unreachable. */
  quadcond: EvidenceResponse | null;
  backend: BackendStatus;
  backendMessage: string;
  selectedId: string;
  overlay: boolean;
  autoRotate: boolean;
  visMode: "ball-and-stick" | "space-filling" | "licorice";
  busy: boolean;
  /**
   * Whether to display values the model itself flagged as out of domain.
   *
   * Default false, and that default is the point. An out-of-domain number is
   * not a slightly-less-certain prediction; it is the model reporting that the
   * query lies outside what it was trained on. Rendered next to an in-domain
   * one, in the same typeface, on the same axis, it is indistinguishable from a
   * real answer. The user can ask to see it -- but has to ask.
   */
  showExtrapolation: boolean;
  setShowExtrapolation: (v: boolean) => void;
  setSequence: (sequence: string) => void;
  setPolymer: (polymer: Polymer) => void;
  setCondition: <K extends keyof Conditions>(key: K, value: Conditions[K]) => void;
  loadExample: (id: string) => void;
  selectStructure: (id: string) => void;
  setOverlay: (v: boolean) => void;
  setAutoRotate: (v: boolean) => void;
  setVisMode: (mode: "ball-and-stick" | "space-filling" | "licorice") => void;
  checkBackend: () => Promise<void>;
  /** Run once on mount: contact the service for the seed sequence. */
  bootstrap: () => void;
}

function toServiceConditions(c: Conditions): QuadcondConditions {
  // A straight copy, now that the UI carries the same eight fields the atlas
  // records. The previous version sent the single "monovalent" slider as
  // potassium and said so in a comment -- which meant a user who typed a
  // sodium buffer got a potassium answer, and the assumption lived in a source
  // file rather than on screen. Across 411 sequence-matched pairs the measured
  // K+/Na+ difference averages +12.9 degC, so that assumption was worth about
  // thirteen degrees of silent error.
  return {
    k: c.k,
    na: c.na,
    li_nh4: c.li_nh4,
    mg: c.mg,
    ph: c.ph,
    temperature: c.temperature,
    crowder_pct: c.crowder_pct,
    strand_conc: c.strand_conc,
  };
}

/**
 * In-flight generation counter.
 *
 * Dragging a slider fires a request per frame and the replies do not
 * necessarily come back in order. Without this, a slow response for pH 5.2
 * can land after the fast one for pH 7.4 and leave the panel showing numbers
 * for a condition the user has already moved away from — which looks exactly
 * like the model being wrong.
 */
/**
 * Map a service card onto a local geometry id.
 *
 * The service answers per structure ("g-quadruplex"); the geometry builder
 * produces one member per drawable topology ("g-quadruplex-parallel",
 * "g-quadruplex-antiparallel"). When the service reports a topology, prefer the
 * matching archetype; otherwise take the first of that kind.
 */
function matchLocalId(qc: EvidenceResponse, result: GeometryResult): string | undefined {
  const card = qc.evidence?.[0];
  if (!card) return undefined;
  return memberForCard(result.members, card)?.id;
}

let generation = 0;
let bootstrapped = false;

export const useNA = create<NAState>((set, get) => {
  const seed = geometryOnly(first.sequence, first.polymer, DEFAULT_CONDITIONS);

  async function refresh(sequence: string, polymer: Polymer, conditions: Conditions) {
    // Local first: geometry and topology parsing are the app's own job and are
    // instant, so the structure redraws immediately while the calibrated
    // numbers are still in flight.
    const local = geometryOnly(sequence, polymer, conditions);
    const mine = ++generation;
    set({
      sequence,
      polymer: local.polymer,
      conditions,
      result: local,
      // Empty, never an invented id. The literals here ("ss-coil",
      // "g-quadruplex") matched no member, and the viewers papered over that
      // with `?? members[0]` — so a selection resolving to nothing quietly drew
      // the first structure instead.
      selectedId: local.members[0]?.id ?? "",
      quadcond: null,
      busy: true,
    });
    try {
      const qc = await fetchEvidence(sequence, toServiceConditions(conditions));
      if (mine !== generation) return; // a newer request already won
      set({
        quadcond: qc,
        backend: "connected",
        backendMessage: `QuadCond ${qc.quadcond_version}, model ${qc.model_version}`,
        // Keep the user's selection pointing at a structure that exists
        // locally. Service card ids ("g-quadruplex") and geometry ids
        // ("g-quadruplex-parallel") are different namespaces, so selecting by
        // the service id silently deselected everything and the overlay drew a
        // different fold from the one the panel highlighted.
        selectedId: matchLocalId(qc, local) ?? local.members[0]?.id ?? "",
        busy: false,
      });
    } catch (err) {
      if (mine !== generation) return;
      set({
        quadcond: null,
        backend: "offline",
        backendMessage:
          err instanceof Error ? err.message : "the QuadCond service is not reachable",
        busy: false,
      });
    }
  }

  return {
    sequence: first.sequence,
    polymer: first.polymer,
    conditions: { ...DEFAULT_CONDITIONS },
    result: seed,
    quadcond: null,
    backend: "unknown",
    backendMessage: "",
    selectedId: seed.members[0]?.id ?? "",
    overlay: false,
    autoRotate: false,
    visMode: "ball-and-stick",
    busy: false,
    showExtrapolation: false,

    setSequence: (sequence) => {
      const { polymer, conditions } = get();
      void refresh(sequence, polymer, conditions);
    },
    setPolymer: (polymer) => {
      const { sequence, conditions } = get();
      void refresh(sequence, polymer, conditions);
    },
    setCondition: (key, value) => {
      const { sequence, polymer, conditions } = get();
      void refresh(sequence, polymer, { ...conditions, [key]: value });
    },
    loadExample: (id) => {
      const ex = EXAMPLES.find((e) => e.id === id);
      if (!ex) return;
      void refresh(ex.sequence, ex.polymer, get().conditions);
    },
    selectStructure: (id) => set({ selectedId: id }),
    setOverlay: (overlay) => set({ overlay }),
    setAutoRotate: (autoRotate) => set({ autoRotate }),
    setVisMode: (visMode) => set({ visMode }),
    setShowExtrapolation: (showExtrapolation) => set({ showExtrapolation }),
    bootstrap: () => {
      // Nothing called the service until the user moved a control, so a fresh
      // page sat on "QuadCond not contacted" with no predictions, indefinitely
      // — and the offline banner claimed the service was unreachable when it
      // had never been asked. Found by rendering the app, which is the only way
      // this class of defect surfaces.
      if (bootstrapped) return;
      bootstrapped = true;
      const { sequence, polymer, conditions } = get();
      void refresh(sequence, polymer, conditions);
    },
    checkBackend: async () => {
      const h = await fetchHealth();
      set({
        backend: h.status === "ok" ? "connected" : "offline",
        backendMessage:
          h.status === "ok" ? `QuadCond ${h.quadcond_version}` : h.error ?? "service unavailable",
      });
    },
  };
});
