"""specs/184: the measured call, the prediction combine and the statistics, on synthetic data."""
import math
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[2]))
import hpn_stim_bench as H  # noqa: E402
import stim_184 as S  # noqa: E402


def line(tmp_path, rows, names=("A_pS1", "B_pS2")):
    p = tmp_path / "X_main.csv"
    head = [",,,Slide ID," + ",".join(names), ",,,Antibody Name," + ",".join(names), ",,,HUGO ID," + ",".join(names),
            "Cell Line,Inhibitor,Stimulus,Timepoint"]
    p.write_text("\n".join(head + [",".join(["X", *map(str, r)]) for r in rows]) + "\n")
    return p


def data(tmp_path, egf_a=4.0, pbs_a=1.0, inhibitor_row=True):
    rows = [("DMSO", "", "0min", 1.0, 1.0), ("DMSO", "", "0min", 1.1, 1.0), ("DMSO", "", "0min", 0.9, 1.0)]
    for t in S.EARLY + S.LATE:
        rows += [("DMSO", "PBS", t, pbs_a, 2.0), ("DMSO", "EGF", t, egf_a, 2.0)]
    if inhibitor_row:
        rows.append(("GSK690693", "EGF", "5min", 100.0, 100.0))
    return S.load(line(tmp_path, rows))


def test_load_keeps_only_dmso_rows_and_takes_log2(tmp_path):
    d = data(tmp_path)
    assert len(d[("", "0min")]) == 3 and len(d[("EGF", "5min")]) == 1
    assert d[("EGF", "5min")][0]["A_pS1"] == pytest.approx(2.0)


def test_measured_is_ligand_over_pbs_at_the_same_time_point(tmp_path):
    m = S.measured(data(tmp_path), "EGF", S.EARLY)
    e, delta = m["A_pS1"]
    assert e == pytest.approx(2.0)
    assert delta > S.FLOOR                    # 2 sd of the time-0 replicates (about 0.29) beats the floor
    assert S.call_of(e, delta) == "up"
    # B is 2.0 under both EGF and PBS: a medium-change effect, not a ligand effect
    assert m["B_pS2"][0] == pytest.approx(0.0) and m["B_pS2"][1] == pytest.approx(S.FLOOR)
    assert S.call_of(*m["B_pS2"]) == "not"


def test_measured_down_and_doubled_threshold(tmp_path):
    d = data(tmp_path, egf_a=0.7)
    assert S.call_of(*S.measured(d, "EGF", S.EARLY)["A_pS1"]) == "down"
    assert S.call_of(*S.measured(d, "EGF", S.EARLY, mult=2.0)["A_pS1"]) == "not"


def test_measured_needs_two_time_points(tmp_path):
    rows = [("DMSO", "", "0min", 1.0, 1.0), ("DMSO", "", "0min", 1.1, 1.0),
            ("DMSO", "PBS", "5min", 1.0, 1.0), ("DMSO", "EGF", "5min", 4.0, 1.0)]
    assert S.measured(S.load(line(tmp_path, rows)), "EGF", S.EARLY) == {}


def test_predictions_combine_networks(tmp_path):
    p = tmp_path / "pred.tsv"
    cols = ["pathway_id", "condition", "readout", "status", "predicted", "reach_predicted", "predicted_withdraw"]
    rows = [("P1", "EGF", "A_pS1", "scored", "up", "up", "not"), ("P2", "EGF", "A_pS1", "scored", "not", "not", "up"),
            ("P1", "EGF", "B_pS2", "scored", "up", "up", "up"), ("P2", "EGF", "B_pS2", "scored", "down", "up", "up"),
            ("P1", "EGF", "C", "trivial", "", "", "")]
    p.write_text("\t".join(cols) + "\n" + "\n".join("\t".join(r) for r in rows) + "\n")
    pred = S.predictions(p)
    assert pred[("EGF", "A_pS1")] == {"model": "up", "reach": "up", "withdraw": "up"}
    assert pred[("EGF", "B_pS2")]["model"] == "not"          # networks disagree: a conflict counts as no change
    assert ("EGF", "C") not in pred


def test_cases_and_unpredicted_count(tmp_path):
    pred = {("EGF", "A_pS1"): {"model": "up", "reach": "up", "withdraw": "not"}}
    cases, unpred = S.build_cases({"X": data(tmp_path)}, pred, ("EGF",), S.EARLY)
    assert [(c["readout"], c["measured"], c["model"]) for c in cases] == [("A_pS1", "up", "up")]
    assert unpred == 1


def case(measured, model, reach, cell="X", ligand="EGF"):
    return {"cell": cell, "ligand": ligand, "readout": "r", "measured": measured, "model": model, "reach": reach,
            "never": "not"}


def test_arm_stats_precision_recall_direction():
    cs = [case("up", "up", "up"), case("up", "down", "up"), case("not", "up", "up"), case("down", "not", "not")]
    s = S.arm_stats(cs, "model")
    assert (s["calls"], s["right"]) == (3, 1)
    assert s["precision"] == pytest.approx(1 / 3) and s["recall"] == pytest.approx(1 / 3)
    assert (s["direction"], s["direction_n"]) == (0.5, 2)
    assert S.arm_stats(cs, "never")["precision"] is None


def test_enrichment():
    cs = [case("up", "up", "up"), case("not", "up", "up"), case("up", "not", "not")] + [case("not", "not", "not")] * 3
    assert S.enrichment(cs) == pytest.approx(0.5 / 0.25)
    assert S.enrichment([case("up", "up", "up")]) is None


def test_bootstrap_resamples_whole_conditions():
    good = [case("down", "down", "up", ligand="EGF")] * 5 + [case("not", "not", "up", ligand="EGF")] * 5
    same = [case("up", "up", "up", ligand="HGF")] * 5
    b = S.bootstrap(good + same, S.ba_diff("model", "reach"), lambda c: c["ligand"], n=300)
    assert b["groups"] == 2 and b["point"] > 0
    assert b["ci"][0] == pytest.approx(0.0) and b["ci"][1] > 0      # a draw of HGF alone has no difference
    assert math.isfinite(b["ci"][1])


def test_engagement_reads_the_receptor(tmp_path):
    rows = [("DMSO", "", "0min", 1.0, 1.0), ("DMSO", "", "0min", 1.1, 1.0)]
    for t in S.EARLY:
        rows += [("DMSO", "PBS", t, 1.0, 1.0), ("DMSO", "EGF", t, 4.0, 1.0), ("DMSO", "HGF", t, 1.0, 1.0)]
    d = S.load(line(tmp_path, rows, names=("EGFR_pY1068", "c-Met_pY1235")))
    e = S.engagement({"X": d}, S.EARLY)
    assert e[("X", "EGF")] is True and e[("X", "HGF")] is False and e[("X", "IGF1")] is None


def test_bench_readout_table_and_helpers():
    ros = {r["readout"]: r for r in H.readouts()}
    assert "PKC-delta" not in ros and len(ros) == 47
    assert ros["EGFR_pY1068"]["sites"] == "EGFR:1092"
    assert ros["4EBP1_pT37_pT46"]["sites"] == "EIF4EBP1:37,EIF4EBP1:46"
    assert H.condition_roots({"EGF": ["b", "a"], "HGF": ["c"], "INS": None}, ["EGF", "INS", "HGF"]) == ["a", "b", "c"]
    assert H.OPPOSITE[H.bm.DOWN] == H.bm.UP and H.OPPOSITE[H.bm.NORMAL] == H.bm.NORMAL
