"""The live-Atlas response parser, against both obs schemas the client emits.

This exists because the first version knew only one of them and answered the
other with an empty dict -- which downstream reads as "no regulatory effect at
this locus". A parser that cannot read a response must say so.
"""
import pytest

from quadcond import alphagenome as ag

anndata = pytest.importorskip("anndata")
np = pytest.importorskip("numpy")
pd = pytest.importorskip("pandas")


class _Variant:
    """What the live client puts in obs.variant."""

    def __init__(self, chrom, pos, ref, alt):
        self.chromosome, self.position = chrom, pos
        self.reference_bases, self.alternate_bases = ref, alt


def _frame(obs, matrix):
    return anndata.AnnData(X=np.asarray(matrix, dtype=float), obs=obs,
                           var=pd.DataFrame(index=[f"t{i}" for i in range(np.shape(matrix)[1])]))


def test_variant_object_schema_is_read_and_aggregated_over_tracks():
    obs = pd.DataFrame({"variant": [_Variant("chr8", 127735930, "G", "T"),
                                    _Variant("chr8", 127735931, "G", "A")]})
    recs = ag._records_from_anndata({"AVI_SCORE": _frame(obs, [[0.4, -1.2], [0.1, 0.05]])}, "live")
    scores = {(k.position, k.alternate): v.scores["AVI_SCORE"] for k, v in recs.items()}
    # the aggregate is max |value| across tracks, so the -1.2 wins as 1.2
    assert scores == {(127735930, "T"): 1.2, (127735931, "A"): 0.1}


def test_four_column_schema_still_works():
    obs = pd.DataFrame({"chromosome": ["chr8"], "position": [127735930],
                        "reference_bases": ["G"], "alternate_bases": ["T"]})
    recs = ag._records_from_anndata({"AVI_SCORE": _frame(obs, [[0.9]])}, "live")
    assert [v.scores for v in recs.values()] == [{"AVI_SCORE": 0.9}]


def test_rows_that_cannot_be_read_raise_rather_than_look_quiet():
    obs = pd.DataFrame({"mystery": [1, 2]})
    with pytest.raises(ag.AtlasUnavailable):
        ag._records_from_anndata({"AVI_SCORE": _frame(obs, [[0.0], [0.0]])}, "live")


def test_no_frames_is_not_an_error():
    assert ag._records_from_anndata({}, "live") == {}
