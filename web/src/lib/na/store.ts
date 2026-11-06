import { create } from "zustand";
import { EXAMPLES } from "./examples";
import { geometryOnly, memberForCard, type GeometryResult } from "./geometry-only";
import {
  DEFAULT_CONDITIONS,
  type Conditions,
  type Polymer,
} from "./types";
import {
  QuadcondRefused,
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
  /**
   * True while `quadcond` belongs to the previous request and a new one is in
   * flight. The panel keeps the last answer on screen, dimmed, instead of
   * blanking to "waiting" on every keystroke and slider step.
   */
  stale: boolean;
  /**
   * The service is running and declined this input (RNA, too short, empty).
   * Not an outage: `backend` stays "connected".
   */
  serviceRefusal: string | null;
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
  /** Several condition fields at once, as one request (e.g. a head's training buffer). */
  applyConditions: (patch: Partial<Conditions>) => void;
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

/**
 * Wait this long after the last change before asking the service. A slider
 * drag used to send one request per pointer move; the generation counter kept
 * the answer right, but every intermediate request still blanked the panel.
 */
export const REQUEST_DEBOUNCE_MS = 200;
let wakePending: (() => void) | null = null;
function settle(ms: number): Promise<void> {
  if (wakePending) wakePending(); // release the superseded wait; its generation check returns
  if (ms <= 0) return Promise.resolve();
  return new Promise((resolve) => {
    const t = setTimeout(() => {
      wakePending = null;
      resolve();
    }, ms);
    wakePending = () => {
      clearTimeout(t);
      wakePending = null;
      resolve();
    };
  });
}

/** Sequence + polymer of the last request the service answered. */
let answeredFor = "";

/**
 * Keep the user's selected structure across a refresh when only the buffer
 * changed. Lowering the pH used to snap the 3D view back to the first card —
 * the G4 — while the user was looking at the i-motif.
 */
export function keepsSelection(
  selectedId: string,
  previousKey: string,
  key: string,
  members: { id: string }[],
): boolean {
  return previousKey === key && !!selectedId && members.some((m) => m.id === selectedId);
}

export const useNA = create<NAState>((set, get) => {
  const seed = geometryOnly(first.sequence, first.polymer, DEFAULT_CONDITIONS);

  async function refresh(sequence: string, polymer: Polymer, conditions: Conditions, delay = REQUEST_DEBOUNCE_MS) {
    // Local first: geometry and topology parsing are the app's own job and are
    // instant, so the structure redraws immediately while the calibrated
    // numbers are still in flight.
    const local = geometryOnly(sequence, polymer, conditions);
    const mine = ++generation;
    const key = `${local.polymer}:${sequence.trim().toUpperCase()}`;
    const prev = get();
    set({
      sequence,
      polymer: local.polymer,
      conditions,
      result: local,
      // Empty, never an invented id. The literals here ("ss-coil",
      // "g-quadruplex") matched no member, and the viewers papered over that
      // with `?? members[0]` — so a selection resolving to nothing quietly drew
      // the first structure instead.
      selectedId: keepsSelection(prev.selectedId, answeredFor, key, local.members)
        ? prev.selectedId
        : (local.members[0]?.id ?? ""),
      // The previous answer stays on screen, marked stale, until this one lands.
      stale: prev.quadcond != null,
      busy: true,
    });
    if (!sequence.trim()) {
      // Nothing to ask about. Not a refusal and not an outage: an empty box.
      answeredFor = key;
      set({ quadcond: null, stale: false, serviceRefusal: null, busy: false });
      return;
    }
    await settle(delay);
    if (mine !== generation) return; // superseded while waiting
    try {
      const qc = await fetchEvidence(sequence, toServiceConditions(conditions));
      if (mine !== generation) return; // a newer request already won
      const current = get().selectedId;
      const keep = keepsSelection(current, answeredFor, key, local.members);
      answeredFor = key;
      set({
        quadcond: qc,
        stale: false,
        serviceRefusal: null,
        backend: "connected",
        backendMessage: `QuadCond ${qc.quadcond_version}, model ${qc.model_version}`,
        // Keep the user's selection pointing at a structure that exists
        // locally. Service card ids ("g-quadruplex") and geometry ids
        // ("g-quadruplex-parallel") are different namespaces, so selecting by
        // the service id silently deselected everything and the overlay drew a
        // different fold from the one the panel highlighted.
        selectedId: keep ? current : (matchLocalId(qc, local) ?? local.members[0]?.id ?? ""),
        busy: false,
      });
    } catch (err) {
      if (mine !== generation) return;
      if (err instanceof QuadcondRefused) {
        // The service answered: there is no prediction for this input. The
        // service is not down, so nothing on screen may say it is.
        answeredFor = key;
        set({
          quadcond: null,
          stale: false,
          serviceRefusal: err.message,
          backend: "connected",
          busy: false,
        });
        return;
      }
      set({
        quadcond: null,
        stale: false,
        serviceRefusal: null,
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
    stale: false,
    serviceRefusal: null,
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
    applyConditions: (patch) => {
      const { sequence, polymer, conditions } = get();
      void refresh(sequence, polymer, { ...conditions, ...patch }, 0);
    },
    loadExample: (id) => {
      const ex = EXAMPLES.find((e) => e.id === id);
      if (!ex) return;
      void refresh(ex.sequence, ex.polymer, get().conditions, 0);
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
      void refresh(sequence, polymer, conditions, 0);
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
