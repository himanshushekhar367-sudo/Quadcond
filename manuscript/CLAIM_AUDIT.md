# Manuscript claim audit

All file references below are relative to the paired release folder. No model fitting was performed in this revision.

| Manuscript location / claim | Evidence anchor | Interpretation retained |
| --- | --- | --- |
| Abstract and Table 1 regression metrics | `quadcond/artifacts/quadcond_model.json`, `heads.*.metrics` | Frozen grouped internal results, not new validation |
| Condition-feature ablation and Figure 2 | Main and `quadcond_model_seqonly.json` sidecars | Same development framework; no error bars or significance claim |
| Atlas size and evidence tiers | `atlas_snapshot.summary` in main sidecar; Table S2 | 398,375 records; experimental tier includes genomic labels |
| Training rows, groups and seeds | `heads.*.training_meta`; Table S1 | Rows are not independent sequences; grouping mitigates but does not prove absence of family leakage |
| Two-construct i-motif condition models | Training metadata, source ingestion and 5DUVMA primary study | C9 and Tel21C only; condition-group holdout is not unseen-sequence validation |
| Empirical uncertainty coverage | `split_conformal_coverage` and `split_conformal_note` | Held-out-group residual partitions; no conditional or mutation-difference guarantee |
| Genomic association interpretation | Service/client labels, Zanin et al., antibody and CUT&Tag literature | Proxy association, not physical folding probability or same-molecule occupancy |
| Refusal rendering | `client.ts`, `EvidencePanels.tsx`, browser results | Refused response types omit numeric values and bypass extrapolation |
| Correct comparison accounting | `quadcond/scans.py`, `tests/test_release_integration.py` | Internally assigned indices and bounded unique mutation positions |
| Readiness configuration | `quadcond/assets.py`, explicit-path regression tests | Readiness and inference inspect the same primary asset configuration |
| Schematic 3D | Geometry code, card matching and independent gate policy | Structural-class illustration; no atomic coordinate prediction |
| Fresh software test counts | `quadcond/artifacts/pytest-results.xml`, frontend verification logs | Software behaviour only; synthetic fixture is not a trained model |
| Public deployment | No verified public endpoint or archive | Explicit placeholders, no claim of availability |
| Independent benchmarking | Existing benchmark artifacts show limitations | No new independent predictive or mutation-effect validation claimed |
| Table 3 scope comparison | Each row restates the question its cited tool describes itself as answering | Scope only; no tool was run, no output was compared, and no ranking is stated |
| Condition dependence of folding (Introduction) | Bhattacharyya et al. on cation identity; Deep et al. on polyamine disruption of i-motifs | Cited as established condition dependence, not as validation of any head |
| Calibration and applicability method citations | Niculescu-Mizil & Caruana; Zadrozny & Elkan; Netzeva et al.; Mitchell et al. | Cited as the standard practice followed, not as evidence that this model is well calibrated outside its training task |

## References

`verified_references.json` retains retrieval URLs and PubMed/Europe PMC metadata, plus Crossref metadata for XGBoost. The manuscript contains 32 numbered references in first-citation order. Every in-text reference resolves; all bibliography entries are cited. DOI links and PMC or PubMed links are included where available. XGBoost is linked to the primary conference DOI. `reference_audit.json` records the exact reference order and, per key, which store the metadata came from.

The 14 references added in this revision were retrieved from the scite metadata store, because Europe PMC was unreachable from the machine that assembled them. `reference_audit.json` marks them and lists the five whose volume or page the retrieved record did not carry. Those locators are left empty rather than guessed: a locator invented to make a reference list look finished is a fabricated locator, and it is the kind of detail nobody re-checks at proof stage.

The unrelated molecular-evolution result returned during one preliminary DOI lookup — Easteal (1985), reached through `10.1093/oxfordjournals.molbev.a040361` — was correctly excluded from the earlier draft. The intended source was Altschul & Erickson's dinucleotide-preserving permutation method, whose DOI is `10.1093/oxfordjournals.molbev.a040370`; it is now cited where the shuffled negatives are described, and the wrong DOI has been corrected in `fetch_references.py` so the mistake is not re-fetched.

References were selected for the actual claims made. Their inclusion is not evidence that comparator software was benchmarked.
