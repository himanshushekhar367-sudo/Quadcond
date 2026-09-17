"""Fetch Atlas scores and run QuadCond on the same GRCh38 motif window.

Run without arguments for prompts. --export-only saves scores without loading
a model. --reuse EXPORT_FOLDER runs QuadCond on an existing export, without a key.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
DEFAULT_ROOT = (HERE.parent if (HERE.parent / 'quadcond' / 'variants.py').is_file()
                else HERE.parent / 'quadcond-repo')


def prompt_value(label, convert, explanation):
    while True:
        value = input(label).strip()
        try:
            return convert(value)
        except ValueError:
            print(explanation)


def chromosome_value(value):
    value = value.strip()
    if value.startswith('chr'):
        value = value[3:]
    if value not in {str(i) for i in range(1, 23)} | {'X', 'Y'}:
        raise ValueError('Use a human chromosome chr1-chr22, chrX or chrY')
    return 'chr' + value


def start_value(value):
    value = int(value)
    if value < 1:
        raise ValueError('Start must be positive')
    return value


def sequence_value(value):
    value = value.strip().upper()
    if not 1 <= len(value) <= 200 or set(value) - set('ACGT'):
        raise ValueError('Enter 1-200 DNA bases using only A, C, G and T')
    return value


def load_quadcond(root):
    if root:
        root = Path(root).resolve()
        if not (root / 'quadcond' / 'variants.py').is_file():
            raise ValueError(f'Not a QuadCond source directory: {root}')
        sys.path.insert(0, str(root))
    import quadcond
    if quadcond.__version__ not in {'0.4.9', '0.5.0', '0.5.1'}:
        raise ValueError(f'Expected QuadCond 0.4.9, 0.5.0 or 0.5.1, found {quadcond.__version__}')
    from quadcond import assets, variants
    return assets, variants


def validate_export(folder):
    """Reject malformed, duplicate or mismatched join keys before model loading."""
    for name in ('provenance.json', 'quadcond_scores.csv'):
        if not (folder / name).is_file():
            raise ValueError(f'No complete Atlas export in {folder}: missing {name}. '
                             'Run --export-only with a chromosome, start and sequence first. '
                             'Use --reuse only after the export reports success.')
    meta = json.loads((folder / 'provenance.json').read_text(encoding='utf-8'))
    if meta.get('assembly') != 'GRCh38':
        raise ValueError('Expected GRCh38 provenance')
    sequence = meta['sequence']
    start = meta['start_1based']
    if not isinstance(start, int) or start < 1 or not sequence or set(sequence) - set('ACGT'):
        raise ValueError('Invalid sequence or start in provenance')
    if meta['interval_0based'] != [start - 1, start - 1 + len(sequence)]:
        raise ValueError('Inconsistent coordinate conventions in provenance')
    import math
    seen = set()
    with (folder / 'quadcond_scores.csv').open(encoding='utf-8-sig', newline='') as f:
        for row in csv.DictReader(f):
            pos = int(row['position'])
            key = (row['chromosome'], pos, row['ref'], row['alt'])
            offset = pos - start
            if key in seen:
                raise ValueError('Duplicate variant in score table')
            seen.add(key)
            if row['chromosome'] != meta['chromosome'] or not 0 <= offset < len(sequence):
                raise ValueError('Score table contains a different genomic window')
            if row['ref'] != sequence[offset] or row['alt'] not in {'A', 'C', 'G', 'T'} or row['alt'] == row['ref']:
                raise ValueError('Score table allele does not match sequence')
            if row['scorer'] != meta['scorer'] or not math.isfinite(float(row['score'])):
                raise ValueError('Mixed scorers or non-finite score in table')
    if not seen:
        raise ValueError('No scored variants in table')
    if len(seen) != meta['returned_scored_snv_count']:
        raise ValueError('Score table row count disagrees with provenance')
    return meta


def model_path(assets, explicit, root):
    if explicit or os.environ.get('QUADCOND_MODEL'):
        return assets.resolve('model', explicit)
    try:
        return assets.resolve('model')
    except assets.AssetError:
        # Also inspect the supplied release and adjacent local artifacts.
        candidates = [Path(root) / 'artifacts' / 'quadcond_model.joblib'] if root else []
        candidates.append(HERE.parent / 'artifacts' / 'quadcond_model.joblib')
        errors = []
        for candidate in candidates:
            if candidate.is_file():
                try:
                    return assets.resolve('model', str(candidate))
                except assets.AssetError as exc:
                    errors.append(str(exc))
        raise assets.AssetError('\n'.join(errors) or 'No matching trained model found; provide --model PATH')


def flatten(result):
    rows = []
    for v in result['variants']:
        s, r = v['structural'], v.get('regulatory') or {}
        rows.append(dict(**v['variant'], structural_head=s['head'],
            structural_delta=s['delta'], structural_strand=s['strand'],
            motif_state=s['motif_state'], regulatory_scorer=r.get('scorer'),
            regulatory_score=r.get('score'), quadrant=v['quadrant'],
            combined_rank=v.get('combined_rank'),
            applicability_json=json.dumps(s.get('applicability'), default=str),
            all_heads_json=json.dumps(v['heads'], default=str)))
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--chromosome')
    p.add_argument('--start', type=int)
    p.add_argument('--sequence')
    p.add_argument('--scorer', default='AVI_SCORE')
    p.add_argument('--ontology', action='append')
    p.add_argument('--out', type=Path, default=HERE / 'motif_export')
    p.add_argument('--reuse', type=Path, help='Existing exporter output directory')
    p.add_argument('--export-only', action='store_true')
    p.add_argument('--quadcond-root', default=str(DEFAULT_ROOT) if DEFAULT_ROOT.exists() else None)
    p.add_argument('--model', help='Trained model; must match release manifest checksum')
    p.add_argument('--head', default='g4_tm', choices=['g4_tm', 'im_pht'])
    p.add_argument('--k', type=float, default=100)
    p.add_argument('--na', type=float, default=0)
    p.add_argument('--mg', type=float, default=0)
    p.add_argument('--ph', type=float, default=7.0)
    p.add_argument('--temperature', type=float, default=25)
    a = p.parse_args()
    if a.reuse and (a.chromosome or a.start or a.sequence or a.ontology):
        p.error('--reuse takes its locus and tissue filters from saved provenance')
    if a.reuse:
        folder = a.reuse.resolve()
    else:
        if not a.chromosome:
            a.chromosome = prompt_value('Chromosome (e.g. chr22): ', chromosome_value,
                'A chromosome is required. Type chr1-chr22, chrX or chrY; press Ctrl+C to cancel.')
        else:
            a.chromosome = chromosome_value(a.chromosome)
        if a.start is None:
            a.start = prompt_value('First base of motif window (GRCh38, 1-based): ', start_value,
                'Type the positive integer coordinate of your sequence\'s FIRST base. This cannot be blank.')
        else:
            a.start = start_value(a.start)
        if not a.sequence:
            a.sequence = prompt_value('Exact forward-strand DNA sequence of motif window: ', sequence_value,
                'Paste 1-200 A/C/G/T bases from the GRCh38 forward strand. This cannot be blank.')
        else:
            a.sequence = sequence_value(a.sequence)
        folder = a.out.resolve()
        command = [sys.executable, str(HERE / 'export_atlas.py'),
            '--chromosome', a.chromosome, '--start', str(a.start),
            '--sequence', a.sequence, '--scorer', a.scorer, '--out', str(folder)]
        for term in a.ontology or []:
            command.extend(['--ontology', term])
        subprocess.run(command, check=True)
    meta = validate_export(folder)
    if a.export_only:
        print(f'QuadCond-compatible export ready: {folder}')
        return
    result_json, result_csv = folder / 'joined_variants.json', folder / 'joined_variants.csv'
    if result_json.exists() or result_csv.exists():
        raise ValueError('Joined outputs already exist. Preserve them and use a new export folder for another run.')
    assets, variants = load_quadcond(a.quadcond_root)
    path = model_path(assets, a.model, a.quadcond_root)
    from quadcond.models.predict import Predictor
    from quadcond.conditions import Condition
    from quadcond.alphagenome import TableAtlas
    from quadcond.service import _check_provenance
    pred = Predictor.load(path, None)
    _check_provenance(pred.model)
    if a.head not in pred.model.heads:
        raise ValueError(f'Model has no {a.head} head')
    condition = Condition.from_mapping(dict(k=a.k, na=a.na, mg=a.mg,
        ph=a.ph, temperature=a.temperature), track_imputed=True)
    print(f'Running {a.head}: K={a.k}, Na={a.na}, Mg={a.mg} mM; pH={a.ph}; T={a.temperature} C')
    result = variants.variant_scan(pred, meta['sequence'], meta['chromosome'],
        meta['start_1based'], condition=condition, heads=[a.head], structural_head=a.head,
        regulatory_scorer=meta['scorer'],
        atlas_source=TableAtlas.from_path(folder / 'quadcond_scores.csv'))
    if result['run_record'].get('regulatory_error'):
        raise RuntimeError(result['run_record']['regulatory_error'])
    result['bridge_provenance'] = dict(atlas=meta, model_path=str(path),
        model_sha256=assets.sha256(path),
        table_sha256=hashlib.sha256((folder / 'quadcond_scores.csv').read_bytes()).hexdigest(),
        interpretation='Two prediction axes; not validated causal or clinical evidence')
    result_json.write_text(json.dumps(result, indent=2, default=str), encoding='utf-8')
    from export_atlas import write_csv
    write_csv(result_csv, flatten(result))
    print(f'Saved joined comparison: {result_csv}')
    print(f'Full predictions, flags and provenance: {result_json}')


if __name__ == '__main__':
    try:
        main()
    except (Exception, KeyboardInterrupt) as exc:
        print(f'Run stopped: {exc}', file=sys.stderr)
        print('Existing files were retained. --reuse requires an export that already completed successfully.', file=sys.stderr)
        sys.exit(1)
