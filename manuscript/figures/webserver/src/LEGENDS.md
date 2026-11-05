# Web-server figures: legends

Numbers in these legends were read from the server's responses when the figures
were made (QuadCond 0.5.1, model file SHA-256 `493640193ad70974…`, which matches
the v0.5.1 release manifest).

## Figure WS1. Architecture of the QuadCond web server

A single process (`quadcond serve`) serves the AENNA-3D workbench and the HTTP API
on one port, so the browser calls the API on its own origin and sends no
credentials. The server has no accounts or sessions and stores no requests.
Requests are checked before scoring: a sequence containing U is refused as RNA,
`/predict` accepts up to 64 sequences of up to 5,000 nt each, and `/scan/genome`
accepts up to 20 kb. The 12 prediction heads are grouped by what their training
label measures: seven are biophysically anchored, three are genomic proxies built
on CUT&Tag peaks, and two learn labels derived from the sequence or predicted by
another model. Each head applies its own applicability gate and returns a
supported estimate, hides an extrapolation, or refuses without a number. Every
response carries a provenance block giving the model file's SHA-256, whether it
matches the release manifest, and how the runtime environment differs from the
one the model was built in. Model and atlas files are downloaded from the GitHub
release and checked by SHA-256 against `assets_manifest.json` when the server
image is built, so a mismatched file stops the build instead of reaching a user.
AlphaGenome variant scores are optional, and are joined to a mutation scan only
on exact chromosome, position, reference and alternate allele.

## Figure WS2. The AENNA-3D interface

**(A)** The workbench with the 21-nt human telomeric repeat
GGGTTAGGGTTAGGGTTAGGG in the default cytosol-like buffer (140 mM K⁺, 10 mM Na⁺,
1 mM Mg²⁺, pH 7.4, 37 °C). Left, sequence entry and buffer controls; centre, a
schematic 3D archetype drawn at published helical parameters (not a predicted
atomic structure) and the Compare tools; right, one evidence card per head.
**(B)** The `g4_tm` card for this strand: predicted Tm 71.9 °C, with a 90%
interval of 59.5–84.3 °C whose half-width is taken from out-of-fold residuals;
the folded fraction of the isolated strand under a two-state model; and the
statement that this is not occupancy at a genomic locus. **(C)** The same strand
with K⁺ removed and Na⁺ raised to 150 mM: predicted Tm 58.7 °C (46.3–71.1 °C).
**(D)** The `im_pht` card for the telomeric C-strand CCCTAACCCTAACCCTAACCCT in
the buffer its training measurements share (100 mM K⁺, 10 mM Na⁺, no Mg²⁺,
pH 7.0, 25 °C): predicted transitional pH 6.38 (5.75–7.01). **(E)** The same
C-strand in the default buffer of (A). Every training record for this head used
100 mM K⁺, no Mg²⁺, pH 7 and 25 °C, so the head has no learned response to any of
the four values requested; it declines, names each one, and withholds the number
unless the user chooses to show the extrapolation. **(F)** The Conditions tool
after sweeping K⁺ from 0 to 150 mM for the strand in (A). `g4_tm` is reported as
responding (its prediction moves 33.5 °C against a typical error of 12.4 °C, a
descriptive comparison rather than a significance test), while heads with no
supported points are listed as such rather than drawn. Screenshots were taken
from a fresh clone of the public repository serving the v0.5.1 release assets.

## Figure WS3. A condition sweep returned by `/scan/conditions`

Predicted Tm (`g4_tm`) of the 21-nt human telomeric repeat as K⁺ (teal, no Na⁺)
or Na⁺ (orange, no K⁺) is swept from 0 to 150 mM in 10 mM steps, with no Mg²⁺,
pH 7.0, 25 °C and the default 5 µM strand concentration. Shaded bands are the 90%
intervals the server reports; their half-width (12.4 °C) is the 90th percentile of
the out-of-fold residuals, so they are not conformal intervals, and the two bands
overlap at every concentration. All 32 points fell inside the head's
applicability domain. Predicted Tm is higher in K⁺ than in Na⁺ at every
concentration above zero, by 7.5–16.8 °C (9.3 °C at 100 mM). Neighbouring points
are not always monotonic (K⁺, 10 → 20 mM) because the ensemble has no
monotonicity constraint, so the endpoints carry the information and the local
slopes do not. The evidence that the head separates the two cations comes from
the held-out cation-identity check (411 K⁺/Na⁺ pairs over 151 sequences; mean
out-of-fold ΔTm +13.6 °C against +12.9 °C measured). This figure shows what the
server returns for one sequence and is not additional validation.
