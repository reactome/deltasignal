import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import score_124  # noqa: E402

COLS = ["pathway", "gene", "direction", "key_output", "predicted", "expected", "valid", "pred_ui", "converged"]


def write(d, axis, rows):
    d.mkdir(parents=True, exist_ok=True)
    with open(d / f"{axis}_cases.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, COLS, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(dict(zip(COLS, r)))


def run(tmp, ctrl_rows, arm_rows, mild_rows, ceiling=1000.0, exp=None):
    for name, rows in (("c", ctrl_rows), ("a", arm_rows), ("m", mild_rows)):
        write(tmp / name, "curator", rows)
        write(tmp / name, "experimental", exp[name] if exp else [])
    return score_124.arm_report(tmp / "c", tmp / "a", ceiling, tmp / "m")


def test_d0_identical_and_different(tmp_path):
    rows = [("Signaling_by_EGFR", "EGFR", "0", "1", "0", "0", "1", "0.200000", "1")]
    for n in ("x", "y"):
        write(tmp_path / n, "curator", rows)
        write(tmp_path / n, "experimental", [])
    out = score_124.d0(tmp_path / "x", tmp_path / "y")
    assert out["curator"] == {"cases": 1, "only_control": 0, "only_frozen": 0}
    write(tmp_path / "y", "curator", [("Signaling_by_EGFR", "EGFR", "0", "1", "0", "0", "1", "0.200001", "1")])
    assert score_124.d0(tmp_path / "x", tmp_path / "y")["curator"]["pred_ui"] == 1


def test_changes_split_by_direction_and_split(tmp_path):
    held, tune = "Signaling_by_EGFR", "Signaling_by_WNT"
    ctrl = [(held, "A", "0", "1", "0", "0", "1", "0.10", "1"),      # right, stays
            (held, "B", "2", "1", "2", "2", "1", "99.5", "1"),      # right at the rail, becomes normal: broke
            (held, "C", "2", "2", "1", "2", "1", "1.00", "0"),      # wrong, becomes up: fixed
            (tune, "D", "0", "1", "0", "1", "1", "0.50", "1")]      # tuning split, unchanged
    arm = [(held, "A", "0", "1", "0", "0", "1", "0.10", "1"),
           (held, "B", "2", "1", "1", "2", "1", "1.10", "1"),
           (held, "C", "2", "2", "2", "2", "1", "3.00", "1"),
           (tune, "D", "0", "1", "0", "1", "1", "0.50", "1")]
    mild = [(held, "B", "2", "1", "2", "2", "1", "99.9", "1"),     # rail-dependent
            (held, "C", "2", "2", "1", "2", "1", "1.00", "1")]
    out = run(tmp_path, ctrl, arm, mild)["curator"]["all"]
    h = out["held-out"]
    assert (h["cases"], h["changed"], h["fixed"], h["broke"]) == (3, 2, 1, 1)
    assert h["changed_rail_dependent"] == 1 and h["changed_mild_known"] == 2
    assert h["control_railed"] == 1 and h["arm_railed"] == 0
    assert h["moved"] == 2 and h["control_converged"] == 2 and h["arm_converged"] == 3
    assert out["held-out|dir2"]["changed"] == 2 and out["held-out|dir0"].get("changed", 0) == 0
    assert out["tuning"]["cases"] == 1 and out["tuning"].get("changed", 0) == 0
    assert h["mcnemar_p"] == 1.0


def test_invalid_cases_and_raf_set_aside(tmp_path):
    ctrl = [("RAF_MAP_kinase_cascade", "BRAF", "0", "1", "2", "0", "1", "50", "1"),
            ("Signaling_by_EGFR", "E", "0", "1", "0", "0", "0", "0.1", "1")]
    arm = [("RAF_MAP_kinase_cascade", "BRAF", "0", "1", "0", "0", "1", "0.2", "1"),
           ("Signaling_by_EGFR", "E", "0", "1", "1", "0", "1", "1.0", "1")]
    rep = run(tmp_path, ctrl, arm, [])["curator"]
    assert rep["all"]["tuning"]["fixed"] == 1 and "held-out" not in rep["all"]   # the invalid case is not scored
    assert "tuning" not in rep["raf_set_aside"]
    assert rep["per_pathway_net"] == {"RAF_MAP_kinase_cascade": 1}


def test_rails_use_each_arms_ceiling_and_zero():
    assert score_124.railed(9.95, 10.0) and not score_124.railed(9.85, 10.0)
    assert score_124.railed(990.0, 1000.0) and not score_124.railed(99.5, 1000.0)
    assert score_124.railed(0.0, 100.0) and score_124.railed(0.01, 100.0) and not score_124.railed(0.02, 100.0)


def test_moved_and_spearman():
    assert not score_124.moved(0.0, 0.0) and score_124.moved(0.0, 0.5)
    assert score_124.moved(2.0, 2.3) and not score_124.moved(2.0, 2.1)
    assert abs(score_124.spearman([1, 2, 3, 4], [10, 20, 30, 40]) - 1.0) < 1e-12
    assert abs(score_124.spearman([1, 2, 3], [3, 2, 1]) + 1.0) < 1e-12
    assert score_124.ranks([5, 5, 1]) == [1.5, 1.5, 0.0]
    assert math.isnan(score_124.spearman([1, 1, 1], [1, 2, 3]))
