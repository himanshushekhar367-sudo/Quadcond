"""Run the supplied chr22 example without coordinate or sequence prompts."""
import datetime
import json
from pathlib import Path
import subprocess
import sys
import requests
import uuid

HERE = Path(__file__).resolve().parent


def fetch_example():
    # User's example variant: GRCh38 chr22:36201698 A>C.
    # Retrieve 30 reference bases on each side, on the forward strand.
    position = 36201698
    start0, end0 = position - 1 - 30, position + 30
    url = ('https://api.genome.ucsc.edu/getData/sequence?'
           f'genome=hg38;chrom=chr22;start={start0};end={end0}')
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    data = response.json()
    sequence = data['dna'].upper()
    if (data.get('genome') != 'hg38' or data.get('chrom') != 'chr22'
            or data.get('start') != start0 or data.get('end') != end0):
        raise ValueError('Reference service returned a different region')
    if len(sequence) != 61 or set(sequence) - set('ACGT'):
        raise ValueError('Reference sequence is incomplete or contains ambiguous bases')
    if sequence[30] != 'A':
        raise ValueError('Reference allele is not A; stopping to avoid a mismatched variant')
    return sequence, start0 + 1, url


def main():
    print('Example: GRCh38 chr22:36201698 A>C', flush=True)
    print('Fetching the surrounding 61 DNA bases automatically...', flush=True)
    sequence, start, url = fetch_example()
    print(f'Reference checked. Window starts at {start}.', flush=True)
    print('This example tests data integration; it is not a confirmed G4/i-motif locus.', flush=True)
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    folder = HERE / 'results' / f'chr22_example_{stamp}_{uuid.uuid4().hex[:6]}'
    subprocess.run([sys.executable, str(HERE / 'run_quadcond.py'), '--export-only',
        '--chromosome', 'chr22', '--start', str(start), '--sequence', sequence,
        '--out', str(folder)], check=True)
    (folder / 'reference_source.json').write_text(json.dumps(dict(
        url=url, assembly='GRCh38/hg38', strand='+', sequence=sequence,
        start_1based=start, example_variant='chr22:36201698:A:C',
        note='Example context only, not an experimentally established structural motif'),
        indent=2), encoding='utf-8')
    (folder / 'reference.fasta').write_text(
        f'>hg38_chr22_{start}_{start + len(sequence) - 1}_forward\n{sequence}\n', encoding='utf-8')
    print('\nSUCCESS: Atlas scores are saved in QuadCond-compatible format.')
    print(f'Open this folder: {folder}')
    print('Main file: quadcond_scores.csv')
    print('Structural model prediction has not been run; the matching model is still required.')


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print('\nCancelled. Existing results were preserved.', file=sys.stderr)
        sys.exit(1)
    except Exception as exc:
        print(f'Could not complete the example: {exc}', file=sys.stderr)
        print('Existing results were preserved. Copy this error into the chat.', file=sys.stderr)
        sys.exit(1)
