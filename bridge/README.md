# Easy AlphaGenome -> QuadCond -> AENNA integration

## QuadCond 0.5.0: use the isolated model environment

From the `quadcond-repo` directory, create this environment once:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e . -r requirements-model.txt
```

Start the workbench with your existing export (no AlphaGenome key or new query):

```powershell
.\.venv\Scripts\python.exe bridge\serve_table.py ..\alphagenome_bridge\results\chr22_example_20260915_105922_ce722e\quadcond_scores.csv
```

Open http://localhost:8765/ and use Compare > Variant join. For that export,
the chromosome is `chr22`, first base `36201668`, and sequence:

```text
AGGAGGGCAGAGAGCTGGGGCCTCGGACTCACCCGACGCTTGTGATGAGCTGCACCCAGGA
```

The export has 183 AVI-scored substitutions. It is a negative-control software
case: missing structural deltas remain missing and variants stay unclassified.
Use the environment's Python for model inference; the machine's default Python
may contain incompatible versions. Stop the old server before starting another
on port 8765.

## New: one-command workflow

Run in PowerShell:

```powershell
cd C:\Users\pc\Downloads\Quadcond\alphagenome_bridge
python run_quadcond.py --export-only
```

Enter your motif's chromosome, GRCh38 first-base coordinate, complete forward
DNA sequence and hidden API key when prompted. This saves a validated
QuadCond-compatible export in `motif_export`. The scorer list you obtained
confirms AVI_SCORE is available for your account.

To fetch and immediately run the full comparison, omit `--export-only`.
This requires the matching trained model. The runner automatically imports
the current `quadcond-repo` source, so a separate editable install is not
required when its Python dependencies are already present.

To run an already exported window without another API request:

```powershell
python run_quadcond.py --reuse motif_export --head g4_tm --k 100 --ph 7 --temperature 25
```

For i-motif comparisons select `--head im_pht`. Conditions are explicit example
values; set them to your experimental conditions. The run records imputed
condition fields. Results are `joined_variants.csv` and `joined_variants.json`
inside the export directory. The JSON retains full predictions and provenance;
the CSV includes structural deltas, Atlas score, applicability, motif state
and the existing QuadCond prioritization labels. Existing joined files are
never overwritten. A partial Atlas table does not classify missing SNVs as low.

Current local model check: the release expects SHA-256
`204b4608c964793107ec5a384581b316c4bd0aacb20791aeb5561df9187a57db`.
The adjacent local model hashes to
`3ada393c2071a3bb6ee11b50524fecf3d1e00f7aa6eda4ba7edf3329d4c16a29`;
the combined sandbox parts hash to
`493640193ad70974d35582c55b4a721095bfcba50edc6827d9636d2f7b7be0c7`.
Neither matches. Supply the matching artifact using `--model PATH`; the
runner preserves the release's checksum check and will stop until it matches.
The supplied manifest has no model download URL. Model execution remains
unverified, and the checksum should not be changed merely to accept a file.

Run software tests with `python -m unittest test_bridge -v`. They use synthetic
test inputs solely to check SDK conversion, reference/coordinate errors,
duplicate/mixed-scorer rejection, the actual QuadCond join (with a mocked
structural predictor), missing-score handling and bad-model rejection.

The v0.4.9 release already contains the variant join, CSV reader and AENNA panel.
Use the Atlas API to retrieve a small subset, save it locally, then reuse it.
No retraining or full Atlas download is required.

## 1. Check access (PowerShell, outside the Python >>> prompt)

```powershell
cd C:\Users\pc\Downloads\Quadcond\alphagenome_bridge
python export_atlas.py --list-scorers
```

The script asks for your API key with hidden input. It does not save the key.
Your installed AlphaGenome version is 0.9.0 and already has the Atlas client.
Use the exact scorer names returned by this command. The exporter defaults to
AVI_SCORE but refuses it if it is not in the service metadata.

An optional one-base access check using the first variant in your transcript:

```powershell
python export_atlas.py --chromosome chr22 --start 36201698 --sequence A --out access_check
```

This retrieves all three alternate alleles at that base, if available. It is
only an access check, not a useful G4/i-motif input window. A reference mismatch
means you must check assembly, coordinate and allele before continuing.

## 2. Export your actual motif window

Choose a short window spanning the complete candidate motif and any desired
flanks. Use exact GRCh38 forward-strand sequence and its first-base coordinate.
The exporter permits up to 200 bases (600 possible SNVs). That computational
limit is not a claim about QuadCond's scientific applicability.

```powershell
$chromosome = Read-Host 'Chromosome, e.g. chr22'
$start = [int](Read-Host 'GRCh38 first base of window, 1-based')
$sequence = Read-Host 'Exact forward-strand DNA sequence of that window'
python export_atlas.py --chromosome $chromosome --start $start --sequence $sequence --out motif_export
```

Outputs:

- `quadcond_scores.csv`: chromosome, position, ref, alt, scorer, score.
- `atlas_details.csv`: signed raw score, available quantile score, all returned
  gene/observation and track metadata. Use this to interpret tissue and direction.
- `provenance.json`: assembly, window, SDK version, query time, scorer, filters,
  aggregation rule and returned/expected SNV counts.

The summary is maximum absolute raw score within the selected scorer over all
returned tracks/genes. For a scalar AVI result this simply retains its magnitude.
It does not combine different scorer types. Missing values remain absent.
For tissue-specific molecular effects, choose a molecular scorer from
`--list-scorers` and optionally add `--ontology YOUR_ONTOLOGY_ID`. Keep these
separate from AVI. The summary alone cannot establish increased/decreased
expression or a causal structural mechanism.

The exporter checks reference bases at returned positions against your pasted
sequence. Missing Atlas positions are not independently sequence-verified.

## 3. Use the existing QuadCond command

If QuadCond 0.4.9 is not installed in this Python environment:

```powershell
python -m pip install -e C:\Users\pc\Downloads\Quadcond\release-v0.4.9\quadcond
```

Use the release's verified model assets and existing asset setup. The inspected
release artifacts directory has model metadata JSON but no trained `.joblib`
file. The sandbox contains split model files from an older package; do not
assume they match the v0.4.9 asset manifest. This adapter does not change model
assets, install dependencies, or bypass checks.

Once your QuadCond model installation is operational:

```powershell
Remove-Item Env:ALPHAGENOME_AVI_TABIX -ErrorAction SilentlyContinue
quadcond variant $sequence --chromosome $chromosome --start $start --atlas-table motif_export\quadcond_scores.csv --regulatory-scorer AVI_SCORE --k 100 --na 0 --mg 0 --ph 7.0 --temperature 25 | Out-File -Encoding utf8 joined_variants.json
```

Those buffer values are an explicit example, not a physiological or validated
default. Choose the actual experimental conditions and a relevant structural
head. Inspect applicability/refusal flags before interpreting any delta. If you
exported another scorer, replace AVI_SCORE with that exact name. The environment
variable removal prevents the original source resolver from preferring a Tabix
file over your specified table; it affects only this PowerShell session.

## 4. Connect to AENNA's existing Variant join panel

Instead of your usual QuadCond backend process, start:

```powershell
python serve_table.py motif_export\quadcond_scores.csv
```

This serves on http://127.0.0.1:8765. Keep it running and connect the existing
AENNA interface to that local backend. Enter the SAME sequence, chromosome and
first-base coordinate in Variant join, choose conditions, and click Join.
The table only covers the exported window; export again for a different locus.
Restart this service to load a different table. Only run one backend on port
8765. The launcher supports the inspected 0.4.9 and 0.5.0 cached-source interface and is
version-guarded. It preserves normal model provenance/asset checks.

A hosted AENNA site may need a reachable backend with its normal deployment
configuration; this launcher is a local prototype and does not deploy a site.

## Why your original call failed

`predict_variant` requires `ontology_terms`; it generates tracks rather than
looking up Atlas scores. Supported context lengths in your installed SDK are
16384, 131072, 524288 and 1048576 bases, so 10000/20000-base intervals are also
unsuitable. The separate `AtlasClient.query_interval` accepts a short query
region and returns precomputed variant scores.

Corrected prediction example, if you want RNA-seq tracks separately:

```python
from alphagenome.models import dna_client
from alphagenome.data import genome

# model = dna_client.create(API_KEY), as in your existing session
variant = genome.Variant('chr22', 36201698, 'A', 'C')
interval = variant.reference_interval.resize(1048576)
result = model.predict_variant(
    interval=interval,
    variant=variant,
    requested_outputs=[dna_client.OutputType.RNA_SEQ],
    ontology_terms=['UBERON:0002107'],  # liver; choose relevant available tissue
)
print(result.reference.rna_seq.values.shape)
print(result.alternate.rna_seq.values.shape)
```

These arrays are predicted reference/alternate tracks, not AVI. The earlier
`atlas_scores.csv` contains error messages, not successful scores. Client
creation alone does not demonstrate a successful authenticated query.

## Verification and scientific limits

Tested locally with actual AlphaGenome Variant objects and AnnData containers:
signed detail preservation, repeated gene aggregation, quantiles, missing
values, sparse matrices, off-by-one/wrong-chromosome/reference rejection, and
CSV round-trip through the exact TableAtlas reader in release-v0.4.9.zip.

Reproduced the release live-reader bug: `_records_from_anndata` silently returns
an empty mapping for the installed SDK's `obs.variant` layout. This adapter
bypasses that reader. No original archive or application code was modified.

Live API authentication/query, trained-model execution and browser end-to-end
behavior have NOT been tested. No API key was supplied in this task. Mock data
were used only for software checks and were not delivered as scientific scores.

Keep the two outputs as separate axes: a predicted structural change and an
AlphaGenome effect. Agreement supports experimental prioritization, not proven
structure-mediated regulation. AENNA remains a schematic viewer. Improved
ranking or publishability needs comparative validation against each tool alone.

Sources checked:
- https://www.alphagenomedocs.com/api/atlas.html
- https://www.alphagenomedocs.com/api/generated/alphagenome.models.dna_client.DnaClient.html
- https://github.com/google-deepmind/alphagenome
- Local AlphaGenome 0.9.0 source and supplied release-v0.4.9.zip.
