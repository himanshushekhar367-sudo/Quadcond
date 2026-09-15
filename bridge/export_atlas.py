"""Small GRCh38 Atlas query -> detailed CSV + QuadCond score table.

Run --help. API key is read from ALPHAGENOME_API_KEY or a hidden prompt.
No model inference, genome download, or QuadCond retraining is required.
"""
import argparse
import csv
import datetime
import getpass
import importlib.metadata
import json
import os
from pathlib import Path

import numpy as np
from alphagenome.atlas import atlas
from alphagenome.data import genome


def convert(frames, chromosome, start, sequence):
    """Preserve track/gene rows; aggregate only within each named scorer.

    The table uses maximum absolute raw score over returned genes/tracks.
    This is an explicitly lossy prioritization summary, never an AVI substitute.
    """
    details, grouped = [], {}
    for scorer, adata in frames.items():
        matrix = adata.X.toarray() if hasattr(adata.X, 'toarray') else np.asarray(adata.X)
        quantiles = adata.layers.get('quantiles')
        if hasattr(quantiles, 'toarray'):
            quantiles = quantiles.toarray()
        if 'variant' not in adata.obs:
            raise ValueError('Unexpected SDK schema: obs.variant is absent')
        for i, (_, obs) in enumerate(adata.obs.iterrows()):
            v = obs['variant']
            offset = v.position - start
            if v.chromosome != chromosome or not 0 <= offset < len(sequence):
                raise ValueError('Atlas returned a variant outside the requested window')
            if sequence[offset] != v.reference_bases:
                raise ValueError(f'Reference mismatch at {v.chromosome}:{v.position}; check GRCh38, strand and start')
            if v.alternate_bases not in 'ACGT' or len(v.alternate_bases) != 1 or v.alternate_bases == v.reference_bases:
                raise ValueError('Expected an SNV')
            key = (v.chromosome, v.position, v.reference_bases, v.alternate_bases, scorer)
            for j in range(matrix.shape[1]):
                value = float(matrix[i, j])
                q = float(quantiles[i, j]) if quantiles is not None else None
                details.append(dict(chromosome=v.chromosome, position=v.position,
                    ref=v.reference_bases, alt=v.alternate_bases, scorer=scorer,
                    raw_score=value if np.isfinite(value) else '',
                    quantile_score=q if q is not None and np.isfinite(q) else '',
                    observation_metadata=json.dumps({k: str(x) for k, x in obs.items() if k != 'variant'}),
                    track_metadata=json.dumps({k: str(x) for k, x in adata.var.iloc[j].items()})))
                if np.isfinite(value):
                    grouped[key] = max(grouped.get(key, 0.0), abs(value))
    summary = [dict(zip(('chromosome', 'position', 'ref', 'alt', 'scorer', 'score'), (*k, v)))
               for k, v in sorted(grouped.items())]
    if not summary:
        raise ValueError('No finite scores returned; this is not evidence of no effect')
    return details, summary


def write_csv(path, rows):
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--chromosome')
    p.add_argument('--start', type=int, help='1-based position of first sequence base, GRCh38')
    p.add_argument('--sequence', help='Exact GRCh38 forward-strand local motif window, 1-200 bases')
    p.add_argument('--scorer', default='AVI_SCORE', help='Exact name from --list-scorers')
    p.add_argument('--ontology', action='append', help='Optional tissue ontology ID; repeat for multiple tissues')
    p.add_argument('--list-scorers', action='store_true')
    p.add_argument('--out', default='atlas_export', help='New output directory; existing directories are refused')
    a = p.parse_args()
    if not a.list_scorers:
        if not a.chromosome or not a.start or not a.sequence:
            p.error('--chromosome, --start and --sequence are required')
        a.sequence = a.sequence.upper().strip()
        if a.start < 1 or not 1 <= len(a.sequence) <= 200 or set(a.sequence) - set('ACGT'):
            p.error('Use a positive 1-based start and 1-200 A/C/G/T bases')
        if not a.chromosome.startswith('chr'):
            a.chromosome = 'chr' + a.chromosome
        if Path(a.out).exists():
            p.error('Choose a new --out directory to avoid stale or overwritten results')
    key = os.environ.get('ALPHAGENOME_API_KEY') or getpass.getpass('AlphaGenome API key (hidden): ')
    client = atlas.create(key, timeout=60)
    metadata = client.scorer_metadata()
    if a.list_scorers:
        print('\n'.join(sorted(metadata)))
        return
    if a.scorer not in metadata:
        raise ValueError(f'Scorer {a.scorer!r} unavailable. Run --list-scorers and choose an exact name.')
    interval = genome.Interval(a.chromosome, a.start - 1, a.start - 1 + len(a.sequence))
    frames = client.query_interval(interval, requested_scorers=[a.scorer],
        ontology_terms=a.ontology, progress_bar=False, max_workers=1)
    details, summary = convert(frames, a.chromosome, a.start, a.sequence)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=False)
    write_csv(out / 'atlas_details.csv', details)
    write_csv(out / 'quadcond_scores.csv', summary)
    provenance = dict(assembly='GRCh38', chromosome=a.chromosome, start_1based=a.start,
        sequence=a.sequence, interval_0based=[interval.start, interval.end],
        scorer=a.scorer, ontology_terms=a.ontology, source='AlphaGenome Atlas API',
        sdk_version=importlib.metadata.version('alphagenome'),
        retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        aggregation='maximum absolute raw score over returned tracks/genes, within this scorer only',
        expected_snv_count=3 * len(a.sequence), returned_scored_snv_count=len(summary),
        missing_scores='absent, never zero-filled', model_release='not supplied by this adapter')
    (out / 'provenance.json').write_text(json.dumps(provenance, indent=2), encoding='utf-8')
    print(f'Saved {len(summary)} scored SNVs out of {3 * len(a.sequence)} possible to {out.resolve()}')


if __name__ == '__main__':
    main()
