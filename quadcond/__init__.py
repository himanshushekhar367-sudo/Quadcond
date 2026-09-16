"""QuadCond -- condition-aware, calibrated prediction for G-quadruplexes and i-motifs.

    from quadcond import Condition, Predictor

    pred = Predictor.load("artifacts/quadcond_model.joblib", "data/atlas.db")
    pred.predict("GGGTTAGGGTTAGGGTTAGGG", Condition.preset("physiological"))

Three ideas hold the package together:

1. **Conditions are inputs, not footnotes.**  Cation identity and
   concentration, pH, temperature and crowding enter the model explicitly, and
   every prediction states whether the requested condition is inside the region
   the model actually learned.
2. **Probabilities must be calibrated.**  Heads are calibrated on out-of-fold
   predictions with cross-fitted evaluation; regression heads carry
   split-conformal intervals.  Ranking metrics alone are not accepted as
   evidence a score can be used as a probability.
3. **Evidence is retrievable.**  Every atlas record keeps its tier
   (experimental / derived / predicted), method, and DOI, and every prediction
   can be accompanied by the nearest real measurements.
"""

__version__ = "0.5.0"

from .conditions import PRESETS, Condition, condition_distance
from .motifs import Element, find_g4, find_im, find_im_graph, g4hunter_mean, scan
from .thermo import folded_fraction_ph, folded_fraction_thermal

__all__ = [
    "__version__",
    "Condition", "PRESETS", "condition_distance",
    "Element", "find_g4", "find_im", "find_im_graph", "g4hunter_mean", "scan",
    "folded_fraction_ph", "folded_fraction_thermal",
]


def __getattr__(name):  # lazy, so importing the package stays cheap
    if name == "Predictor":
        from .models.predict import Predictor
        return Predictor
    if name == "Atlas":
        from .atlas import Atlas
        return Atlas
    raise AttributeError(name)
