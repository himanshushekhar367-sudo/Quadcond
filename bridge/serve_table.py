"""Attach an exported score table to the QuadCond 0.4.9 / 0.5.x AENNA service.

This small version-specific adapter uses the service's existing cached source
slot. It does not modify the installed package or disable its asset checks.
"""
import argparse
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('table', type=Path)
    p.add_argument('--port', type=int, default=8765)
    a = p.parse_args()
    import quadcond
    from quadcond import service
    from quadcond.alphagenome import TableAtlas
    if not (quadcond.__version__ == '0.4.9' or quadcond.__version__.startswith('0.5.')):
        raise RuntimeError('This launcher supports QuadCond 0.4.9 and 0.5.x only')
    source = TableAtlas.from_path(a.table)
    if not source._records:
        raise RuntimeError('The table has no variants')
    # Set once, before any request threads start. The API key is not needed.
    service._regulatory = source
    print(f'Atlas table: {a.table.resolve()}')
    print('Use the same GRCh38 sequence and first-base coordinate in AENNA Variant join.')
    service.serve(host='127.0.0.1', port=a.port)


if __name__ == '__main__':
    main()
