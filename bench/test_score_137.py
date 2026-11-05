"""specs/137: the matched-cutoff rescoring."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import score_137 as s  # noqa: E402


def row(pw, ds, ro, ui, pred, meas, sec=None):
    return (pw, ds, ro), {"pathway_id": pw, "drug_set": ds, "readout": ro, "status": "scored", "pred_ui": str(ui),
                          "predicted": pred, "measured": meas, "measured_secondary": sec or meas}


def test_call_at_is_symmetric_in_log_and_inclusive():
    assert s.call_at(1.5, 1.5) == "up"
    assert s.call_at(1 / 1.5, 1.5) == "down"
    assert s.call_at(0.8, 1.5) == "not"
    assert s.call_at(1.2, 1.5) == "not"
    assert s.call_at(0.8, 1.15) == "down"
    assert s.call_at(0, 1.5) == "down"     # a predicted zero is below every cutoff


def test_rescored_changes_only_predicted_of_scored_rows():
    k, r = row("P", "d1", "A", 0.8, "down", "not")
    k2, r2 = row("P", "d1", "B", 0.2, "down", "down")
    r3 = dict(r2, status="trivial", readout="C", pred_ui="", predicted="")
    cases = {k: r, k2: r2, ("P", "d1", "C"): r3}
    out = s.rescored(cases, 1.5)
    assert out[k]["predicted"] == "not" and out[k2]["predicted"] == "down"
    assert out[("P", "d1", "C")] == r3
    assert cases[k]["predicted"] == "down"   # input untouched


def test_auroc_perfect_reversed_and_ties():
    assert s.auroc([3, 4], [1, 2]) == 1.0
    assert s.auroc([1, 2], [3, 4]) == 0.0
    assert s.auroc([1], [1]) == 0.5
    assert s.auroc([], [1]) is None


def test_report_gate_and_raf_set_aside():
    cases = dict([row("P", "d1", "A", 0.8, "down", "not"),     # fixed by t = 1.5
                  row("P", "d2", "A", 0.3, "down", "down"),    # stays right
                  row("P", "d3", "A", 1.0, "not", "not"),
                  row(s.RAF, "d1", "A", 0.8, "down", "down"),  # broken by t = 1.5, RAF
                  row("P", "d4", "B", 2.0, "up", "up")])
    out = s.report(cases, 1.5)
    g = out["all"]["gate"]
    assert g["fixed"] == 1 and g["broken"] == 1
    assert out["raf_set_aside"]["gate"]["fixed"] == 1 and out["raf_set_aside"]["gate"]["broken"] == 0
    assert out["raf_set_aside"]["gate"]["control"]["n"] == 4


def test_cv_choice_picks_the_cutoff_the_rest_supports():
    # every drug set: one modest fall measured not, one deep fall measured down -> any t >= 1.25 is perfect
    cases = {}
    for d in ("d1", "d2", "d3"):
        cases.update([row("P", d, "A", 0.8, "down", "not"), row("P", d, "B", 0.1, "down", "down")])
    cv = s.cv_choice(sorted(cases), cases)
    assert cv["chosen_t_by_cases"] == {"1.25": 6}
    assert cv["cv_balanced"] == 1.0 and cv["cv_gain"] > 0
