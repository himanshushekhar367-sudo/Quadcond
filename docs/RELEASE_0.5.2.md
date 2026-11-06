# Release 0.5.2 — what recording the demo found

0.5.2 changes no model and no asset: it runs on the v0.5.1 release assets,
verified against the same SHA-256s. Every fix below came from driving the
workbench end to end for a recorded demonstration, which is the only way most of
them surface.

## 1. The i-motif card disappeared below pH 7

All 718 transitional-pH records in the atlas store `ph = 7.0`, and every one of
them lists `ph` in `condition_imputed`: a titration sweeps pH, so no single
value describes it and 7.0 is a placeholder. The applicability gate read that
placeholder as "this head was only trained at pH 7" and withheld the i-motif
card at every acidic pH — exactly where i-motifs fold.

The fix is the one the gate already applied to temperature on melting-
temperature heads: **a head's own target is not one of its inputs.** pH is now
exempt for `ph_t` heads (`quadcond/models/predict.py`), and
`capabilities.axis_response` reports it as `predicts`. The requested pH is the
evaluation point of the Hill titration, whose cooperativity (n = 2.5) is stated
on the card as an assumption. The head itself is evaluated at the pH its
training rows carry, so the predicted pH_T can no longer drift with the pH the
user asks about (it read 6.37 above pH 7 and 6.38 below; it is now 6.38 at
every pH). Every other axis is still gated. On the telomeric
C-strand in the head's buffer, the folded fraction now reads 0.03 at pH 7.0,
0.47 at 6.4 and 0.99 at 5.5, in domain.

## 2. A withheld card names the buffer that would bring it back

`applicability.suggested_condition` (prediction schema 1.2.0, an added field)
lists only the axes that must move and to what: a fixed axis to its one training
value, a varied axis clamped into its range. The workbench turns it into one
button — "Use the i-motif model's buffer: K⁺ 100 mM · Mg²⁺ 0 mM · 25 °C" — where
the demo previously needed five manual edits.

## 3. A refused input looked like an outage

An RNA example, an empty box or a too-short sequence returns HTTP 400 from a
healthy service, and the client treated every non-200 the same way: header
"Model unavailable", panel "QuadCond unreachable — geometry only", footer
"Prediction service not connected". Refusals (400/413/422) now raise
`QuadcondRefused`; the status stays "Model connected" and the panel shows "No
prediction for this input" with the service's reason. Before the first answer
the status reads "Connecting…" rather than "unavailable".

## 4. The panel blanked on every keystroke

Each change cleared the cards to "Waiting for QuadCond predictions…" and a
slider drag sent one request per pointer move. Requests now wait 200 ms for the
input to settle, and the previous answer stays on screen, dimmed and marked
`aria-busy`, until the new one lands.

## 5. Smaller fixes

- **Number boxes accept typing.** "25" into a field whose minimum is 4 failed at
  the "2"; the box now keeps the draft, commits valid values, flags invalid ones
  and tidies up on blur or Enter.
- **The 3D selection survives a buffer change.** Lowering the pH no longer snaps
  the view from the i-motif back to the first card.
- **The chr22 variant-join example works without configuration.** When no
  regulatory source is configured, the service serves a bundled table holding
  the AlphaGenome scores for exactly that window's 183 substitutions, labelled
  `example_table`; every other window is reported unscored. A table export can
  now also be configured with `ALPHAGENOME_TABLE`.
- **The 3D view renders on demand.** A static scene was redrawn continuously,
  including behind the Compare tab; with software rendering that was ~190% CPU
  at idle, now ~2%. Auto-rotate still runs a frame loop.
- **The mutation scan says how many substitutions destroy the motif** above the
  table. On a C-rich strand every C sits in a C-tract, so 36 of 66 substitutions
  lose the i-motif and have no delta; read without that line, the small loop
  effects in the table looked like the model ignoring the C-tracts.

## Tests

`tests/test_ph_gate.py` and `tests/test_example_source.py` (backend) and
`web/src/lib/na/fixes.test.ts` (workbench) pin each change.
