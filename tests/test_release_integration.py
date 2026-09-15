"""Regression checks for input identity, bounded work and asset configuration."""
import pytest
from quadcond import assets, scans, service
from quadcond.models.predict import Predictor
from test_comparison_contracts import _sentinel_model, TEL22


def test_batch_assigns_unique_indices_even_when_caller_reuses_them():
    result = scans.batch_predict(Predictor(_sentinel_model()), [
        {"id": "first", "sequence": TEL22, "index": 7},
        {"id": "second", "sequence": TEL22, "index": 7},
        {"id": "empty", "sequence": "", "index": 7},
    ])
    record = result["run_record"]
    assert (record["n_sequences_in"], record["n_sequences_scored"],
            record["n_sequences_excluded"], record["n_rows"]) == (3, 2, 1, 2)
    assert {r["record_index"] for r in result["results"]} == {0, 1}
    assert result["excluded"][0]["record_index"] == 2


@pytest.mark.parametrize("positions", [[0, 0], [0] * 1000, [1.5], [True], [-1]])
def test_invalid_mutation_positions_fail_before_prediction(positions, monkeypatch):
    pred = Predictor(_sentinel_model())
    def forbidden(*args, **kwargs):
        pytest.fail("invalid positions reached prediction")
    pred.predict = forbidden
    with pytest.raises(ValueError):
        scans.mutation_scan(pred, TEL22, positions=positions)
    monkeypatch.setattr(service, "predictor", lambda: pred)
    with pytest.raises(ValueError):
        service.mutation_scan_payload({"sequence": TEL22, "positions": positions})


@pytest.mark.parametrize("kind,env", [("model", "QUADCOND_MODEL"), ("atlas", "QUADCOND_DB")])
def test_status_and_resolution_use_the_same_explicit_configuration(tmp_path, monkeypatch, kind, env):
    path = tmp_path / 'explicit.bin'
    path.write_bytes(b'verified fixture bytes')
    asset = assets.Asset('primary', 'default.bin', assets.sha256(path), path.stat().st_size, kind, 'test')
    optional = assets.Asset('optional', 'other.bin', 'unused', 1, kind, 'test', required=False)
    monkeypatch.setattr(assets, 'assets', lambda: {'primary': asset, 'optional': optional})
    monkeypatch.setattr(assets, 'search_paths', lambda a: [tmp_path / a.filename])
    monkeypatch.setenv(env, str(path))
    assert assets.resolve(kind) == path
    rows = assets.status()
    assert rows[0]['state'] == 'ok' and rows[0]['path'] == str(path)
    assert rows[1]['state'] == 'missing'
    path.write_bytes(b'wrong bytes')
    with pytest.raises(assets.AssetError):
        assets.resolve(kind)
    assert assets.status()[0]['state'] != 'ok'
    monkeypatch.setenv(env, str(tmp_path / 'absent.bin'))
    assert assets.status()[0]['state'] == 'missing'
    with pytest.raises(assets.AssetError):
        assets.resolve(kind)
