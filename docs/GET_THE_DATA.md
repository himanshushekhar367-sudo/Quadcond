# Getting the datasets

**Status: all three biophysical datasets are in.** Every measurement-grounded head now
trains on real data. What follows records how, and what is left.

| dataset | status | what it gave |
|---|---|---|
| G4STAB Supplementary Table 1 | **ingested** | 2,367 records, 2,274 usable Tm across 261 buffers |
| iM-Seeker Supplementary Data Set | **ingested** | 160 measured pH_T values, 147 sequences |
| 5DUVMA (gkag110 SI) | **ingested** | 937 records: 379 Tm + 558 pH_T across 10 ionic strengths |
| GSE220882 iMab CUT&Tag | downloaded, not ingested | needs peak calling + hg38 |
| gkad626 SI (Zanin et al.) | held out | 14 CD-tested oligos, kept as an independent check |

None of these can be fetched from the cloud sandbox: the publisher CDNs sit behind a
Cloudflare bot check and the agent proxy blocks them. Each needs one manual download. A
pre-filled adapter ships for each in `adapters/`.

---

## 1. G4STAB Supplementary Table 1 — 2,382 experimental G4 Tm ✅ ingested

Sequence, Tm, pH, and [K⁺] / [Na⁺] / [Li⁺,NH₄⁺] for 2,382 literature-curated melting
measurements. This is the dataset the condition-aware Tm head was designed around, and the
only one that makes the condition ablation meaningful for stability.

**Article** <https://academic.oup.com/bioinformatics/article/41/10/btaf545/8266696>
**DOI** [10.1093/bioinformatics/btaf545](https://doi.org/10.1093/bioinformatics/btaf545)

1. Open the article, clear the Cloudflare check.
2. Scroll to **Supplementary data** near the bottom (or use the anchor
   `…/8266696#supplementary-data`).
3. Download **`btaf545_supplementary_data.zip`**.
4. Unzip it and find the sheet holding Table 1 — the one with a sequence column, a Tm
   column, and salt columns.

```bash
quadcond atlas ingest g4stab-supp "Dataset (G4STAB) Supplementary Table 1.csv"
```

A dedicated Python adapter handles this one, because the table needed more than column
mapping:

- **Buffers are prose.** `Comments` holds 419 distinct strings; `quadcond.buffers` resolves
  99.5% of them into cation concentrations. 11 rows whose buffer could not be resolved are
  dropped rather than defaulted.
- **93 melting values are not numbers.** 46 are censored bounds (`>90`, `<30`, `low`) and
  43 are biphasic pairs (`51/71`). Both are ingested as folding evidence with `tm` NULL and
  a QC flag, never coerced — a fabricated midpoint in the regression target is worse than a
  missing one.
- **pH 0 means unknown**, and is recovered from the buffer text where that states one.

**Then:** `python scripts/02_train.py && python scripts/03_model_card.py && quadcond report`

Result: `g4_tm` at R² 0.670 (grouped CV, 400 sequence clusters), and the condition ablation
worth **+0.416 R²** — the project's central claim, on measurements.

---

## 2. iM-Seeker Supplementary Data Set — i-motif transitional pH ✅ ingested

120–171 i-motif sequences with measured pH_T in 10 mM sodium cacodylate + 100 mM KCl over
pH 4–8. S1 is literature-derived, S2 is the authors' in-house measurements. This creates
`im_pht` — the first genuinely pH-aware output — and lets the placeholder
`im_architecture` head be retired.

**Paper** *Prediction of DNA i-motifs via machine learning*, NAR 2024, 52(5):2188–2202
**DOI** [10.1093/nar/gkae092](https://doi.org/10.1093/nar/gkae092)
**Supplementary** <https://academic.oup.com/nar/article/52/5/2188/7607873#supplementary-data>
**Authors' figshare** <https://figshare.com/s/e4e72e2e8ceaa0a4fbd6> — this is the deposit
the iM-Seeker README points at for the model pickles; check it for the tables too, since
it needs no institutional access.
**Code** <https://github.com/YANGB1/iM-Seeker>

The file to get is **`gkae092_supplemental_files.zip`** from the NAR supplementary; inside
it, `Supplementary Data Set.xlsx`. Sheet 1 holds 171 constructs with a measured pH_T;
Sheet 2 re-expresses 120 of them as extracted motifs.

```bash
quadcond atlas ingest imseeker "Supplementary Data Set.xlsx"
```

Only Sheet 1 is ingested as measurements — ingesting Sheet 2 as well would enter each pH_T
twice under two spellings of the same sequence and put near-identical rows on both sides of
a split. Sheet 2's detected motif is attached to its parent record as a flag instead. 11
sequences written with the authors' Δ deletion notation are rejected, not stripped.

**One trap worth knowing about.** pH_T is the pH at which the structure is half folded, so
it is tempting to store it as the record's condition pH. Doing that lets the pH_T head read
its own target out of an input feature; it reported R² 0.993 before this was caught, and
0.592 afterwards. A titration does not happen *at* one pH. The condition is the salt, the
label is the midpoint, and there is a regression test for it.

---

## 3. 5DUVMA — condition-resolved i-motif stability ✅ ingested

Tel21C ([C₃TAA]₃C₃) and C9 ([C₉T₃]₃C₉) across pH 2.0–9.0 in 0.15-unit steps at ten ionic
strengths (0.01–1.0 M KCl, Britton–Robinson buffer), 5–95 °C at 0.2 °C/min: 270 + 180
complete denaturation/renaturation curves with T½, pH½ and both transition widths.

**Paper** *High-throughput measurement and prediction of the i-motif DNA stability
landscape*, NAR 2026, 54(4):gkag110
**DOI** [10.1093/nar/gkag110](https://doi.org/10.1093/nar/gkag110)
**Article** <https://academic.oup.com/nar/article/54/4/gkag110/8474393>
**PDF** <https://academic.oup.com/nar/article-pdf/54/4/gkag110/66859574/gkag110.pdf>

The data is in the supporting-information PDF (`gkag110_supplemental_file.pdf`), as
tables rather than a spreadsheet, so the adapter parses the PDF text directly — no manual
transcription, and re-running it reproduces the same 937 records.

```bash
quadcond atlas ingest duvma gkag110_supplemental_file.pdf
```

Six grids are ingested: Tel21C T½ (16 pH × 10 ionic strengths) and pH½ (16 temperatures ×
10), and C9 melting and annealing for both quantities. Cells marked `a` in the source
("quality not sufficient for accurate determination") are skipped, never zero-filled. The
stated ionic strength is entered as K⁺, since the Britton–Robinson buffer is brought to
ionic strength with KCl — an assumption recorded on every row rather than buried.

**Why both #2 and #3, and what happened when they were pooled.** They are mirror images:
iM-Seeker is *sequence breadth* at one buffer, 5DUVMA is *condition depth* for two
sequences. The instinct is to pool them into one pH_T head. That produced **R² −0.257** —
worse than the mean — because grouping by sequence then means holding out an entire
condition surface for a construct the model has never seen.

They are two questions, and the fix was to train two heads: `im_pht` (sequence → pH_T,
grouped by sequence, R² 0.592) and `im_pht_condition` (buffer → pH_T for known constructs,
**grouped by buffer** so whole ionic strengths are held out, R² 0.835). Splitting the
question is what made both answerable.

---

## 4. GSE220882 — iMab / BG4 CUT&Tag ✅ ingested

`GSE220882_RAW.tar` holds **17 bigWigs and no peak calls**. Two things about that decided
the pipeline, and both contradict the obvious plan:

- **MACS2 cannot run on it.** `-f BAMPE` needs alignments; a bigWig has none. That route
  means going back to SRA.
- **SEACR is the right caller and could not be run either** — bash + R + bedtools, none
  available where the data lives. What runs instead is a Python reimplementation of
  SEACR's no-control logic, cut at the peak counts the authors publish in their Table S1.
  The implied top-fraction lands at 0.94–1.04% across all twelve samples, so one threshold
  explains the whole series.

No lift-over is needed: the deposited tracks are already hg38.

```bash
python scripts/06_gse220882_peaks.py \
    --bigwig-dir GSE220882_RAW \
    --genome-star /path/to/STAR_index \
    --blacklist data/genomic/hg38-blacklist.v2.bed \
    --out data/external/downloads/gse220882_peaks.jsonl
```

Sequence comes out of a STAR index rather than a genome download (UCSC and Ensembl are
unreachable here), and the decoding is verified against the index's own splice junctions
before anything is written. Full detail, including the measured genomic background rates
and why the negatives are shuffled *windows* rather than shuffled motifs, is in
[`GSE220882.md`](GSE220882.md).

Ingested as `label_class="genomic_proxy"`, `evidence_tier="experimental"`, HEK293T only by
default — the WDLPS libraries are ~10× shallower and their peaks sit barely above genomic
background. These rows train `im_fold_genomic` and `g4_fold_genomic`; `g4_fold` and
`im_fold` stay pinned to biophysical labels so their meaning does not shift underneath
them.

---

## Handing files back to this session

Attach the downloaded file to the chat, or connect the folder you saved it in through the
desktop app, and the ingest + retrain + report cycle runs from there. Roughly 25 minutes
of compute, most of it the distilled-Tm head.
