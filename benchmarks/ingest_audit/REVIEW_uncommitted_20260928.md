# Review: uncommitted RNA changes found in the working tree, 2026-09-28

Two changes were present in the working tree, timestamped 10:08 and 10:09 UTC,
and were **not** produced by the audit work (10:40–10:41). Neither is committed.
This note reviews them; it adopts neither. No test pins either change, and
adding one would itself be adoption.

| Path | State | Wired in? |
| --- | --- | --- |
| `quadcond/atlas/db.py` | modified, uncommitted | yes — `Atlas.query` signature |
| `quadcond/atlas/ingest/gse296171_rg4.py` | new, untracked | **no** — absent from `ADAPTERS` |

The adapter is inert: `quadcond/atlas/ingest/__init__.py` does not import it, so
nothing can ingest through it by accident.

## 1. `db.py` — a DNA-only default on `Atlas.query`

Adds `nucleic_acids: Sequence[str] | None = ("DNA",)` and filters on
`COALESCE(nucleic_acid, 'DNA')`.

**The reasoning is sound and the placement is right.** `motifs.clean()` maps U to
T, so an RNA insert whose sequence also exists as DNA is indistinguishable
downstream. Without the fence it would surface from `neighbours()` — and so in
every served prediction's evidence block — as measured DNA evidence at
similarity 1.0. All ten call sites reach the atlas through `query()`, including
`_ensure_vectors` (`db.py:311`), which is what `neighbours()` is built on, and
the two paths that matter most:

- `quadcond/models/train.py:581`
- `quadcond/evaluate.py:38`

`COALESCE` correctly reads a NULL `nucleic_acid` as DNA, so rows predating the
column's default are unaffected. Against the current all-DNA atlas the change is
a no-op, which is why the suite does not move.

**Four paths bypass it.** Each is raw SQL over `records` and will include RNA
rows while `query()` excludes them:

| Path | Effect once RNA exists |
| --- | --- |
| `db.py:151` `count()` | `count()` and `len(query())` disagree |
| `db.py:220` `composition()` | composition reports rows training cannot see |
| `db.py:267` `summary()` | same |
| `assets.py:146` atlas fingerprint | see below |

The reporting divergence is arguably correct — composition *should* show
everything — but it is currently silent, and a reader comparing `composition()`
against a trained head's row count will be misled. Whichever way it is resolved,
it should be stated rather than left to be discovered.

**The fingerprint is the one real defect.** `assets.py` hashes an explicit column
list that does **not** include `nucleic_acid`. Two rows identical in every listed
column but differing in molecule hash identically, so the atlas identity does not
distinguish a DNA atlas from one carrying RNA rows. That undermines the
manifest's purpose the moment RNA lands, and it is worth fixing whether or not
the fence is adopted.

**Risk of adopting.** It is a silent behaviour change on a default argument:
every existing caller changes meaning, and a caller that genuinely wants all rows
now gets DNA with no warning. Acceptable given the fence's purpose, but it needs a
test, and the test is what makes it survive the next refactor.

## 2. `gse296171_rg4.py` — the RNA adapter

Not registered, so this is a design review only.

**It is more careful than this audit was on conditions, and it caught a real
error here.** It computes in-reaction concentrations from the pipetting volumes
rather than quoting the stock:

```
monovalent = 375 * 4/22 + 100 * 6/22 = 95.45 mM
Mg2+       =  15 * 4/22              =  2.73 mM   (total, not free)
```

`gse296171_overlap.py` had recorded the nominal 100 mM / 15 mM — overstating
Mg²⁺ by **5.5×** on the one axis a condition-aware model exists to represent.
That constant has been corrected and a test now pins it
(`test_conditions_are_in_reaction_not_stock`). The adapter's figures are
themselves provisional: they derive from the preprint methods, and the adapter
says so. They should be checked against the published methods before either file
is trusted.

It also handles well: the authors' own −2 cut-off, attributed to them rather
than invented; the insert-versus-90-nt-oligo distinction, so nobody assumes a
feature computed on the insert saw the whole molecule; and the pH caveat that
Tris at 40 °C sits ~0.4 units below the stock reading.

**It dissolves three of the four schema blockers this audit reported, by
declining rather than extending.** PDS arms are dropped (no ligand field needed),
alternative normalisations are dropped (no triple-counting), RBNS is dropped
(binding, not folding). That leaves only `rt_stop_score`, and the adapter takes
`folded` via the authors' threshold with the continuous score in `qc_flags`.

**So the `schema_readiness: BLOCKED` gate is over-strict, and that is this
audit's error, not the adapter's.** Two of the reported requirements were
conflated:

- `label_class`. The adapter files rows as `biophysical` with
  `measurement=enzymatic_proxy_rt_stop`. Re-reading `atlas/schema.py`, that is
  correct by the schema's own definition — `biophysical` means a structure
  measured in a defined buffer, which RT-stop is, as against `genomic_proxy` for
  in-cell assays. The concern that motivated a new class was *applicability
  domain*: pH 8.3 and 2.7 mM Mg²⁺ is a region no other atlas row occupies. That
  is a fence for the claims module to enforce, not a label class.
- `ligand` / `ligand_conc`. Only required if PDS rows are ingested. They are not.

`rt_stop_score` remains a genuine enhancement — thresholding at −2 discards the
gradation that makes the dataset valuable — but it is not a blocker.

**The leakage gate is unaffected and still blocks.** The fence controls what
`query()` returns; it does not establish which sequences entered which
evaluation fold, and no fold manifest is persisted. The 77 benchmark collisions
should be read as *potential* exposure conditional on the training path, not as
77 confirmed leaks — with the fence adopted and all training routed through
`query()`, DNA heads cannot see the RNA rows at all.

## Recommendation

1. Fix the `assets.py` fingerprint to include `nucleic_acid`. Independent of
   everything else.
2. Decide the `count()` / `composition()` / `summary()` divergence explicitly and
   document it.
3. If the fence is adopted, commit it **with** a test asserting no training path
   can reach RNA rows.
4. Relax this audit's `schema_readiness` gate to require only `rt_stop_score`,
   and reclassify the condition-region concern as applicability rather than
   label class.
5. Check the in-reaction condition arithmetic against the published methods.
6. Leave the adapter unregistered until the fold manifest exists.
