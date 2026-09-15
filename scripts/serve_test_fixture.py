"""Explicit synthetic fixture for browser integration tests. No fitting or artifact writes.

This runs the production HTTP handler and prediction/scan code with deterministic
sentinel estimators. It cannot satisfy the pinned real-model release identity.
"""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'tests'))
from test_comparison_contracts import _sentinel_model, _head
from quadcond import service
from quadcond.models.predict import Predictor


class ClassSentinel:
    def __init__(self, values):
        self.values = values
    def predict_proba(self, X):
        return np.tile(self.values, (len(X), 1))


model = _sentinel_model()
for head in model.heads.values():
    head.applicability['k'] = {'min': 10, 'max': 100, 'n_unique': 12}
for name, kind in [('im_fold_genomic', 'iM'), ('g4_fold_genomic', 'G4')]:
    model.heads[name] = _head(name, kind, 'binary', 'folded', {'length': {'min': 10, 'max': 201}})
for name, kind, target, classes, values in [
    ('g4_topology', 'G4', 'topology', ['parallel', 'antiparallel', 'hybrid'], [.1, .2, .7]),
    ('locus_peak_overlap_state', 'locus', 'topology', ['neither', 'G4', 'iM', 'both'], [.1, .2, .3, .4]),
]:
    head = _head(name, kind, 'multiclass', target, {'length': {'min': 10, 'max': 201}}, classes=classes)
    head.estimators = [ClassSentinel(values)] * 4
    model.heads[name] = head
model.atlas_snapshot = {'provenance': 'synthetic', 'provenance_note':
    'Deterministic control-flow fixture. No values represent measurements or fitted-model predictions.'}
service.ALLOW_SYNTHETIC_MODEL = True
service._predictor = Predictor(model)
service._predictor.synthetic_banner = service._check_provenance(model)
service.serve('127.0.0.1', 8765)
