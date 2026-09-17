"""Exercise the local table launcher without opening a port or querying an API."""
import importlib.util
import sys
from pathlib import Path

import pytest

import quadcond
from quadcond import service


def _load(name):
    path = Path(__file__).resolve().parents[1] / 'bridge' / f'{name}.py'
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_table_launcher_accepts_current_release_and_retains_real_source(tmp_path, monkeypatch):
    table = tmp_path / 'scores.csv'
    table.write_text('chromosome,position,ref,alt,AVI\nchr22,101,A,C,0.5\n')
    seen = {}
    monkeypatch.setattr(quadcond, '__version__', '0.5.0')
    monkeypatch.setattr(service, '_regulatory', None)
    monkeypatch.setattr(service, 'serve', lambda **kw: seen.update(kw))
    monkeypatch.setattr(sys, 'argv', ['serve_table.py', str(table)])
    _load('serve_table').main()
    assert seen == {'host': '127.0.0.1', 'port': 8765}
    assert service._regulatory.describe()['source'] == 'table'
    assert len(service._regulatory._records) == 1


def test_table_launcher_refuses_unreviewed_version(monkeypatch):
    monkeypatch.setattr(quadcond, '__version__', '99.0.0')
    monkeypatch.setattr(sys, 'argv', ['serve_table.py', 'missing.csv'])
    with pytest.raises(RuntimeError, match='supports QuadCond'):
        _load('serve_table').main()


def test_bridge_uses_current_checkout_instead_of_old_release():
    bridge = _load('run_quadcond')
    assert bridge.DEFAULT_ROOT == Path(__file__).resolve().parents[1]
    assets, _ = bridge.load_quadcond(bridge.DEFAULT_ROOT)
    assert assets.MANIFEST_PATH.resolve().is_relative_to(bridge.DEFAULT_ROOT)
