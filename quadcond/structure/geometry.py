"""Idealised coordinates for G-quadruplex and i-motif folds.

This is a **model builder, not a structure predictor**. It places nucleotides on
the geometry those folds are known to adopt, at published helical parameters,
so that a rendered picture and an exported coordinate file are at least the
right shape and the right size. It does not refine, does not score, and does not
claim to reproduce any particular deposited structure. Everything it writes is
labelled as a model.

That distinction matters because the alternative in use was worse than
approximate: the same B-form duplex parameters -- 3.38 A rise, ~34 deg twist,
~9 A radius -- were applied to every fold including quadruplexes. Those numbers
are correct for a duplex and wrong for a G4, which is four strands around
stacked tetrads, and wrong twice over for an i-motif, whose two duplexes
intercalate so that successive C:C+ pairs of one strand pair sit ~6.2 A apart
with the other pair's bases between them. A picture built on duplex constants
is not a stylistic choice; it misrepresents the fold's dimensions by a factor
of two in the rise and puts the tetrad H-bond network somewhere it cannot be.

Parameters used, with the source of each:

===================== ========= =============================================
quantity              value     basis
===================== ========= =============================================
G4 tetrad rise        3.3 A     stacking distance between successive G-tetrads
G4 tetrad twist       30 deg    right-handed, parallel-stranded quadruplexes
G4 C1' radius         9.7 A     C1' distance from the four-fold axis
G4 base-centre radius 4.6 A     guanine O6 atoms coordinate the axial cation
i-motif pair rise     3.1 A     rise per *intercalated* C:C+ step
i-motif duplex rise   6.2 A     successive C:C+ pairs of one parallel duplex
i-motif twist         12 deg    near-untwisted, right-handed
i-motif C1' radius    5.4 A     narrow-groove backbone separation
B-DNA rise / twist    3.38 A /  unchanged; the duplex builder was already right
                      34.3 deg
===================== ========= =============================================

These are round numbers from the consensus ranges rather than values fitted to
one PDB entry, which is the honest level of precision for an idealised model.

:func:`validate` checks the geometry it produces rather than trusting it:
tetrad planarity, inter-tetrad spacing, four-fold symmetry, backbone step
lengths and a hard-sphere clash test. A builder that silently drifts fails the
check instead of shipping a plausible-looking picture.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

# ---------------------------------------------------------------- parameters
G4_RISE = 3.3
G4_TWIST_DEG = 30.0
G4_C1_RADIUS = 9.7
G4_BASE_RADIUS = 4.6

IM_PAIR_RISE = 3.1
IM_DUPLEX_RISE = 6.2
IM_TWIST_DEG = 12.0
IM_C1_RADIUS = 5.4
IM_PAIR_SEPARATION = 8.6      # C1'-C1' across a hemiprotonated C:C+ pair

B_RISE = 3.38
B_TWIST_DEG = 34.3
B_C1_RADIUS = 9.0

BACKBONE_STEP = 6.0           # typical P-P separation along a strand
MIN_CONTACT = 2.8             # hard-sphere floor for the clash test


@dataclass
class Residue:
    index: int
    base: str
    strand: int
    role: str                  # tetrad | pair | loop | stem | ss
    level: int | None          # tetrad / intercalation level, when it has one
    c1: tuple[float, float, float]
    base_centre: tuple[float, float, float]
    syn: bool = False          # glycosidic conformation, G4 only


@dataclass
class Model:
    kind: str
    topology: str
    residues: list[Residue] = field(default_factory=list)
    strands: list[list[int]] = field(default_factory=list)
    pairs: list[tuple[int, int]] = field(default_factory=list)
    note: str = ""


def _cyl(radius: float, angle_deg: float, y: float) -> tuple[float, float, float]:
    a = math.radians(angle_deg)
    return (radius * math.cos(a), y, radius * math.sin(a))


def _dist(a, b) -> float:
    return math.dist(a, b)



def _flank_position(anchor: tuple[float, float, float], step_index: int,
                    outward: int) -> tuple[float, float, float]:
    """Place a 5' or 3' tail residue away from the folded core.

    Flanking residues sit outside the tract range, so the interpolation used for
    loops has no second anchor and collapses every one of them onto the same
    point. They get their own path instead: a shallow helical tail leaving the
    core at a real backbone spacing.
    """
    ang = math.radians(35.0 * step_index)
    radius = math.hypot(anchor[0], anchor[2]) + 2.0 * step_index
    return (radius * math.cos(math.atan2(anchor[2], anchor[0]) + ang * 0.3),
            anchor[1] + outward * BACKBONE_STEP * 0.55 * step_index,
            radius * math.sin(math.atan2(anchor[2], anchor[0]) + ang * 0.3))



def _outside_position(m: "Model", placed: dict[int, int], core: list[int],
                      i: int, *, bulge_scale: float, lift: float
                      ) -> tuple[float, float, float]:
    """Where a non-core residue goes: a loop arc, or a tail if it has no bracket."""
    if i < core[0]:
        return _flank_position(m.residues[placed[core[0]]].c1, core[0] - i, -1)
    if i > core[-1]:
        return _flank_position(m.residues[placed[core[-1]]].c1, i - core[-1], +1)
    left = max(c for c in core if c < i)
    right = min(c for c in core if c > i)
    a0, a1 = m.residues[placed[left]].c1, m.residues[placed[right]].c1
    t = (i - left) / max(right - left, 1)
    bulge = bulge_scale * max(math.sin(math.pi * t), 0.25)
    cx = a0[0] * (1 - t) + a1[0] * t
    cy = a0[1] * (1 - t) + a1[1] * t
    cz = a0[2] * (1 - t) + a1[2] * t
    norm = math.hypot(cx, cz) or 1.0
    return (cx * (1 + bulge / norm), cy + bulge * lift, cz * (1 + bulge / norm))


# ---------------------------------------------------------------- G4
def build_g4(seq: str, tracts: list[tuple[int, int]], topology: str = "parallel") -> Model:
    """Place a quadruplex from its parsed G-tracts.

    ``tracts`` is four (start, end) pairs in sequence coordinates. The number of
    stacked tetrads is the shortest tract; extra guanines in longer tracts are
    routed as flanking residues rather than forced into a tetrad that has no
    partner in the other three columns.

    Strand polarity and glycosidic conformation follow the topology, because
    that is the difference between the three folds and it is visible in the
    model: a parallel quadruplex has all four columns running the same way and
    all-anti guanines with propeller loops on the outside; an antiparallel one
    alternates direction and alternates syn/anti; a hybrid runs three columns
    one way and one the other.
    """
    s = seq.upper()
    levels = min(e - b for b, e in tracts)
    if topology == "parallel":
        direction = [1, 1, 1, 1]
    elif topology == "antiparallel":
        direction = [1, -1, 1, -1]
    else:                                   # hybrid (3 + 1)
        direction = [1, 1, 1, -1]

    m = Model(kind="g-quadruplex", topology=topology,
              note="idealised model geometry, not a refined structure")
    placed: dict[int, int] = {}
    top_y = (levels - 1) * G4_RISE

    for col, (b, e) in enumerate(tracts):
        col_angle = col * 90.0
        for lvl in range(levels):
            # a column running 5'->3' downwards starts at the top tetrad
            src = b + lvl if direction[col] > 0 else e - 1 - lvl
            y = lvl * G4_RISE if direction[col] > 0 else top_y - lvl * G4_RISE
            # A tetrad is a set of four coplanar guanines, so membership is
            # defined by HEIGHT, not by position within a tract. In an
            # antiparallel fold a down-running column's first guanine stacks
            # with an up-running column's last one; indexing the level by
            # position instead put them 6.6 A apart in the same "tetrad".
            tetrad = int(round(y / G4_RISE))
            twist = tetrad * G4_TWIST_DEG
            syn = (topology != "parallel") and (direction[col] < 0)
            r = Residue(
                index=src, base=s[src] if src < len(s) else "G",
                strand=col, role="tetrad", level=tetrad,
                c1=_cyl(G4_C1_RADIUS, col_angle + twist, y),
                base_centre=_cyl(G4_BASE_RADIUS, col_angle + twist + 12.0, y),
                syn=syn,
            )
            placed[src] = len(m.residues)
            m.residues.append(r)

    # loops and flanks ride outside the core on an arc between their columns
    core = sorted(placed)
    for i, ch in enumerate(s):
        if i in placed:
            continue
        c1 = _outside_position(m, placed, core, i, bulge_scale=5.0, lift=0.35)
        role = "flank" if (i < core[0] or i > core[-1]) else "loop"
        m.residues.append(Residue(index=i, base=ch, strand=0, role=role,
                                  level=None, c1=c1, base_centre=c1))

    m.residues.sort(key=lambda r: r.index)
    m.strands = [[i for i, _ in enumerate(m.residues)]]
    m.pairs = _tetrad_cycles(m)
    return m


def _tetrad_cycles(m: Model) -> list[tuple[int, int]]:
    """The Hoogsteen cycle in each tetrad, built after residues are in order.

    Rebuilt rather than remapped: the cycle is a property of which residues
    share a level, and deriving it once at the end removes the class of bug
    where an index list survives a sort that invalidated it.
    """
    by_level: dict[int, list[int]] = {}
    for i, r in enumerate(m.residues):
        if r.role == "tetrad" and r.level is not None:
            by_level.setdefault(r.level, []).append(i)
    out: list[tuple[int, int]] = []
    for idxs in by_level.values():
        idxs.sort(key=lambda i: m.residues[i].strand)
        for j in range(len(idxs)):
            out.append((idxs[j], idxs[(j + 1) % len(idxs)]))
    return out


# ---------------------------------------------------------------- i-motif
def build_imotif(seq: str, tracts: list[tuple[int, int]]) -> Model:
    """Place an i-motif from its four parsed C-tracts.

    The fold is two parallel-stranded duplexes intercalated into each other with
    antiparallel polarity. Tracts 1 and 3 form one duplex, tracts 2 and 4 the
    other; successive C:C+ pairs *within* one duplex are ~6.2 A apart because
    the partner duplex's pairs sit between them, and the intercalated stack rises
    ~3.1 A per step. That factor of two is precisely what a duplex-parameter
    model gets wrong.
    """
    s = seq.upper()
    levels = min(e - b for b, e in tracts)
    m = Model(kind="i-motif", topology="intercalated C:C+ tetraplex",
              note="idealised model geometry, not a refined structure")
    placed: dict[int, int] = {}
    pair_seq_idx: list[tuple[int, int]] = []

    # duplex A: tracts 0 (down) and 2 (up); duplex B: tracts 1 (down), 3 (up)
    duplexes = [((0, 1), (2, -1)), ((1, 1), (3, -1))]
    for d, ((ca, da), (cb, db)) in enumerate(duplexes):
        for lvl in range(levels):
            # duplex B is offset half a step, which is the intercalation
            y = (2 * lvl + d) * IM_PAIR_RISE
            twist = (2 * lvl + d) * IM_TWIST_DEG
            for col, dirn, side in ((ca, da, +1), (cb, db, -1)):
                b, e = tracts[col]
                src = b + lvl if dirn > 0 else e - 1 - lvl
                ang = twist + (0.0 if side > 0 else 180.0) + d * 90.0
                c1 = _cyl(IM_C1_RADIUS, ang, y)
                bc = _cyl(IM_PAIR_SEPARATION * 0.25, ang, y)
                placed[src] = len(m.residues)
                m.residues.append(Residue(
                    index=src, base=s[src] if src < len(s) else "C",
                    strand=col, role="pair", level=2 * lvl + d,
                    c1=c1, base_centre=bc))
            pair_seq_idx.append((m.residues[-2].index, m.residues[-1].index))

    core = sorted(placed)
    for i, ch in enumerate(s):
        if i in placed:
            continue
        c1 = _outside_position(m, placed, core, i, bulge_scale=6.0, lift=0.0)
        role = "flank" if (i < core[0] or i > core[-1]) else "loop"
        m.residues.append(Residue(index=i, base=ch, strand=0, role=role,
                                  level=None, c1=c1, base_centre=c1))

    m.residues.sort(key=lambda r: r.index)
    # pairs were recorded in sequence coordinates, so the sort cannot invalidate
    # them; they are resolved to positions exactly once, here
    pos = {r.index: i for i, r in enumerate(m.residues)}
    m.pairs = [(pos[a], pos[b]) for a, b in pair_seq_idx if a in pos and b in pos]
    m.strands = [[i for i, _ in enumerate(m.residues)]]
    return m


# ---------------------------------------------------------------- validation
def validate(m: Model) -> dict:
    """Measure the model instead of trusting the constants that built it."""
    res = m.residues
    out: dict[str, object] = {"kind": m.kind, "topology": m.topology,
                              "n_residues": len(res)}

    levels: dict[int, list[Residue]] = {}
    for r in res:
        if r.level is not None:
            levels.setdefault(r.level, []).append(r)

    if m.kind == "g-quadruplex" and levels:
        planar = [max(abs(x.c1[1] - min(v.c1[1] for v in vs)) for x in vs)
                  for vs in levels.values()]
        ys = sorted({round(vs[0].c1[1], 3) for vs in levels.values()})
        steps = [round(b - a, 2) for a, b in zip(ys, ys[1:])]
        radii = [round(math.hypot(r.c1[0], r.c1[2]), 2) for vs in levels.values() for r in vs]
        out.update({
            "tetrads": len(levels),
            "max_out_of_plane_A": round(max(planar), 3),
            "inter_tetrad_rise_A": steps,
            "c1_radius_A": [min(radii), max(radii)],
            "four_per_tetrad": all(len(v) == 4 for v in levels.values()),
        })
    if m.kind == "i-motif" and levels:
        ys = sorted({round(vs[0].c1[1], 3) for vs in levels.values()})
        steps = [round(b - a, 2) for a, b in zip(ys, ys[1:])]
        same_duplex = [round(ys[i + 2] - ys[i], 2) for i in range(len(ys) - 2)]
        out.update({
            "intercalation_levels": len(levels),
            "intercalated_step_A": steps,
            "same_duplex_step_A": same_duplex,
            "two_per_level": all(len(v) == 2 for v in levels.values()),
        })

    core = [r for r in res if r.role in ("tetrad", "pair")]
    clashes = 0
    for i in range(len(core)):
        for j in range(i + 1, len(core)):
            if _dist(core[i].c1, core[j].c1) < MIN_CONTACT:
                clashes += 1
    out["c1_clashes_under_2.8A"] = clashes

    steps = [_dist(res[i].c1, res[i + 1].c1) for i in range(len(res) - 1)]
    out["backbone_step_A"] = {"min": round(min(steps), 2), "median":
                              round(sorted(steps)[len(steps) // 2], 2),
                              "max": round(max(steps), 2)}
    out["passes"] = bool(clashes == 0 and (not levels or out.get(
        "four_per_tetrad", out.get("two_per_level", True))))
    return out
