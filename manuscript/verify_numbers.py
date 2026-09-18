#!/usr/bin/env python3
"""Check every quantitative claim in the rendered manuscript against its source.

Each entry is (substring that must appear in the manuscript, how to recompute it).
A claim that cannot be recomputed from inputs/ does not belong in the manuscript.
"""
import csv, json, re
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TXT = (ROOT / 'QuadCond_AENNA_NAR_manuscript.md').read_text(encoding='utf-8')
BEN = {(r['task'], r['tool'], r['metric']): r
       for r in csv.DictReader((ROOT / 'inputs' / 'benchmark_results.csv').open(encoding='utf-8'))}
SCN = {(r['split'], r['negatives'], r['tool']): float(r['auroc'])
       for r in csv.DictReader((ROOT / 'inputs' / 'scanner_benchmark.csv').open(encoding='utf-8'))}
GW = json.loads((ROOT / 'inputs' / 'genomewide_summary.json').read_text(encoding='utf-8'))
SENS = json.loads((ROOT / 'inputs' / 'genomewide_sensitivity.json').read_text(encoding='utf-8'))
POS = list(csv.DictReader((ROOT / 'inputs' / 'posctrl_summary.csv').open(encoding='utf-8')))


def bench(task_sub, tool, metric):
    for (t, to, m), v in BEN.items():
        if task_sub in t and to == tool and m == metric:
            return float(v['value'])
    raise KeyError(f'{task_sub} | {tool} | {metric}')


CHECKS = [
    ('0.672', bench('G4 melting temperature', 'QuadCond g4_tm (grouped OOF)', 'R2'), 3),
    ('5.56 °C', bench('G4 melting temperature', 'QuadCond g4_tm (grouped OOF)', 'MAE (°C)'), 2),
    ('ρ 0.819', bench('G4 melting temperature', 'QuadCond g4_tm (grouped OOF)', 'Spearman'), 3),
    ('R² 0.478', bench('G4 melting temperature', 'G4STAB architecture retrained (grouped 5-fold)', 'R2'), 3),
    ('7.42 °C', bench('G4 melting temperature', 'G4STAB architecture retrained (grouped 5-fold)', 'MAE (°C)'), 2),
    ('R² 0.826', bench('G4 melting temperature', 'QuadCond g4_tm (random 5-fold, for comparison with published G4STAB protocol)', 'R2'), 3),
    ('0.896', bench('Within-sequence buffer response', 'QuadCond g4_tm (grouped OOF)', 'median per-sequence Spearman'), 3),
    ('ρ 0.537', bench('Single-substitution effect', 'QuadCond g4_tm (Δ of grouped OOF predictions)', 'Spearman'), 3),
    ('80.9%', bench('Single-substitution effect', 'QuadCond g4_tm (Δ of grouped OOF predictions)', 'sign accuracy |ΔTm|≥2°C') * 100, 1),
    ('AUROC of 0.851', bench('i-motif folding', 'QuadCond im_fold (grouped OOF)', 'AUROC'), 3),
    ('0.620 for iM-Seeker', bench('i-motif folding', 'iM-Seeker (folding probability)', 'AUROC'), 3),
    ('AUROC of 0.957', bench('G4 folding vs dinucleotide shuffles', 'QuadCond g4_fold (grouped OOF)', 'AUROC'), 3),
    ('AUROC of 0.942 against random', SCN[('test_heldout_chr', 'random', 'QuadCond genome-scan (G4-seq)')], 3),
    ('0.928 against sequence-matched', SCN[('test_heldout_chr', 'pq', 'QuadCond genome-scan (G4-seq)')], 3),
    ('G4Hunter reached 0.672', SCN[('test_heldout_chr', 'pq', 'G4Hunter (max |25-nt window|)')], 3),
    ('scanner reached 0.928 and 0.775', SCN[('test_mouse', 'pq', 'QuadCond genome-scan (G4-seq)')], 3),
    ('2,707,107', GW['n_snv'], 0),
    ('37.1% of 1,195 loci', GW['kinds']['G4']['within_locus_lost_vs_flank']['frac_motifs_lost_higher'] * 100, 1),
    ('36.2% of 2,254 loci', GW['kinds']['iM']['within_locus_lost_vs_flank']['frac_motifs_lost_higher'] * 100, 1),
    ('gives 37.2% and 37.3%', GW['kinds']['G4']['within_locus_motif_vs_flank']['frac_motifs_higher'] * 100, 1),
    ('−0.59 percentile points', GW['kinds']['G4']['adjusted_rank_model']['terms']['motif_lost']['delta_pctile'] * 100, 2),
    ('−0.97 points', GW['kinds']['iM']['adjusted_rank_model']['terms']['motif_lost']['delta_pctile'] * 100, 2),
    ('+4.42 percentile points', GW['kinds']['G4']['adjusted_rank_model']['terms']['control']['delta_pctile'] * 100, 2),
    ('+9.67 percentile points', SENS['kinds']['G4']['cpg_in_flank']['adjusted']['delta_pctile'] * 100, 2),
    ('+10.11 for i-motif', SENS['kinds']['iM']['cpg_in_flank']['adjusted']['delta_pctile'] * 100, 2),
    ('0.018 for G4', GW['kinds']['G4']['dose_response_spearman_minus_delta_vs_abs_avi']['rho'], 3),
    ('0.013 for i-motifs', GW['kinds']['iM']['dose_response_spearman_minus_delta_vs_abs_avi']['rho'], 3),
    ('at most 1.43 points', max(abs(x) for x in GW['kinds']['G4']['adjusted_rank_model']['terms']['motif_lost']['ci']) * 100, 2),
    ('at most 1.54', max(abs(x) for x in GW['kinds']['iM']['adjusted_rank_model']['terms']['motif_lost']['ci']) * 100, 2),
    ('0.69% and 0.65%', SENS['kinds']['G4']['variance_explained_r2']['cpg_plus_subst'] * 100, 2),
    ('explains 0.068% and 0.022%', SENS['kinds']['G4']['variance_explained_r2']['motif_class'] * 100, 3),
    ('9.38 °C at VEGFA Pu22', float(next(r for r in POS if r['control'] == 'VEGFA Pu22')['max_abs_delta_Tm']), 2),
    ('5.58 °C at KIT c-kit1', float(next(r for r in POS if r['control'] == 'KIT c-kit1')['max_abs_delta_Tm']), 2),
    ('462 substitutions were evaluated', sum(int(r['n_with_delta']) + int(r['n_motif_lost']) for r in POS if r['control'] != 'chr22 negative-control window'), 0),
    ('of which 135 abolished', sum(int(r['n_motif_lost']) for r in POS), 0),
    ('183 substitutions', int(next(r for r in POS if r['control'].startswith('chr22'))['n_substitutions']), 0),
]

fails = []
for claim, value, dp in CHECKS:
    if claim not in TXT:
        fails.append(f'NOT IN TEXT: {claim!r}')
        continue
    nums = [float(x.replace('−', '-').replace(',', '')) for x in re.findall(r'-?−?[\d,]+\.?\d*', claim) if any(c.isdigit() for c in x)]
    if not nums:
        continue
    # Round half up, matching the manuscript builder. Python's round() and
    # '%.Nf' both inherit binary rounding, which turns 0.4775 into 0.477 and
    # would flag correctly-rounded prose as a mismatch.
    q = Decimal(1).scaleb(-dp)
    want = float(Decimal(repr(float(value))).quantize(q, rounding=ROUND_HALF_UP))
    if not any(abs(n - want) < 10 ** -dp / 2 + 1e-9 for n in nums):
        fails.append(f'MISMATCH: {claim!r} -> source says {want}')

print(f'{len(CHECKS)} claims checked, {len(fails)} problems')
for f in fails:
    print(' ', f)
raise SystemExit(1 if fails else 0)
