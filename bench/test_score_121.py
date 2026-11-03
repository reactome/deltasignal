import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import score_121 as s  # noqa: E402


def test_pearson_edge_cases():
    assert s.pearson([1, 2], [1, 2]) is None
    assert s.pearson([1, 1, 1], [1, 2, 3]) == 0.0
    assert s.pearson([1, 2, 3], [2, 4, 6]) == pytest.approx(1.0)


def test_predictions_average_pathways(tmp_path):
    f = tmp_path / "c.tsv"
    f.write_text("pathway_id\tdrug\treadout\tstatus\tpred_ui\tpredicted\treach_predicted\n"
                 "P1\tD\tAkt\tscored\t4\tup\tdown\nP2\tD\tAkt\tscored\t1\tnot\tnot\n"
                 "P3\tD\tAkt\ttrivial\t\t\t\n")
    p = s.predictions([f])
    assert p["fold"][("D", "Akt")] == pytest.approx(1.0)       # mean of log2 4 and log2 1
    assert p["call"][("D", "Akt")] == 0.5 and p["reach"][("D", "Akt")] == -0.5


def test_per_drug_means_samples():
    preds = {"fold": {("D", "a"): 1.0, ("D", "b"): 0.0, ("D", "c"): -1.0}}
    samples = [{"drug": "D", "values": {"a": 1, "b": 0, "c": -1}}, {"drug": "D", "values": {"a": -1, "b": 0, "c": 1}}]
    assert s.per_drug(samples, preds, {"a", "b", "c"})["fold"]["D"] == pytest.approx(0.0)


def test_per_drug_gated_uses_the_samples_cell_line():
    gated = {"fold": {("D", "C1", "a"): 1.0, ("D", "C1", "b"): 0.0, ("D", "C1", "c"): -1.0}}
    samples = [{"drug": "D", "cell_line": "C1", "values": {"a": 2, "b": 0, "c": -2}},
               {"drug": "D", "cell_line": "C2", "values": {"a": -2, "b": 0, "c": 2}}]        # C2 has no gated prediction
    assert s.per_drug_gated(samples, gated, {"a", "b", "c"})["fold"]["D"] == pytest.approx(1.0)
