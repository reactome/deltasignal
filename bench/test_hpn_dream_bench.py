"""specs/112: the HPN-DREAM scorer's own arithmetic, on fixtures (no challenge data)."""
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hpn_dream_bench as h  # noqa: E402


def test_auroc_ties_count_half():
    assert h.auroc([1, 0, 0, 0], [1, 0, 0, 1]) == pytest.approx(0.75)
    assert h.auroc([0, 0, 0], [1, 0, 1]) == pytest.approx(0.5)      # all tied
    assert h.auroc([3, 2, 1], [1, 1, 0]) == pytest.approx(1.0)
    assert h.auroc([1, 2, 3], [1, 1, 0]) == pytest.approx(0.0)
    assert h.auroc([1, 2], [1, 1]) is None                           # no negatives


def test_auroc_matches_the_challenge_roc():
    # DREAMTools adds a ROC point only where the score changes and integrates by trapezoid
    scores = [0.9, 0.5, 0.5, 0.5, 0.1, 0.0]
    classes = [1, 0, 1, 1, 0, 0]
    pts, tp, fp, prev = [(0.0, 0.0)], 0, 0, None
    for s, c in sorted(zip(scores, classes), key=lambda x: -x[0]):
        if prev is not None and s != prev:
            pts.append((fp / 3, tp / 3))
        tp += c; fp += 1 - c; prev = s
    pts.append((fp / 3, tp / 3))
    trapz = sum((x2 - x1) * (y1 + y2) / 2 for (x1, y1), (x2, y2) in zip(pts, pts[1:]))
    assert h.auroc(scores, classes) == pytest.approx(trapz)


def test_species_order_reads_dv_columns(tmp_path):
    d = tmp_path / "x_experimental" / "experimental" / "MIDAS"
    d.mkdir(parents=True)
    (d / "MD_MCF7_main.csv").write_text("TR:MCF7:CellLine,TR:EGF:Stimuli,DA:ALL,DV:AKT_pS473,DV:TAZ_pS89\n1,2,3,4,5\n")
    assert h.species_order(tmp_path, "MCF7") == ["AKT_pS473", "TAZ_pS89"]


def test_site_map_is_complete_and_parses():
    m = h.load_site_map()
    assert len(m) == 48
    assert m["EGFR_pY1068"] == [("EGFR", "Y1092")]       # mature -> canonical numbering
    assert m["PDK1_pS241"] == [("PDPK1", "S241")]
    assert m["PKC-delta"] == []                           # total protein: no site


def test_score_aligns_truth_to_header_and_drops_nan(tmp_path):
    d = tmp_path / "x_experimental" / "experimental" / "MIDAS"
    d.mkdir(parents=True)
    t = tmp_path / "x_SC1A_TrueDescVectors" / "TrueDescVectors"
    t.mkdir(parents=True)
    out = tmp_path / "out"
    out.mkdir()
    sp = ["A", "B", "C", "TAZ"]
    for cell in h.CELL_LINES:
        (d / f"MD_{cell}_main.csv").write_text(",".join(["TR:x"] + [f"DV:{s}" for s in sp]) + "\n")
        for stim in h.STIMULI:
            (t / f"TrueVec_{cell}_{stim}.csv").write_text("1,0,0,NaN\n")
    with open(out / "hpn_predictions.tsv", "w") as f:
        f.write("stimulus\tspecies\tds_score\treach_score\tbest_pathway\n")
        for stim in h.STIMULI:
            for s, v in zip(sp, [2.0, 1.0, 1e-14, 9.0]):
                f.write(f"{stim}\t{s}\t{v}\t{int(v >= 1)}\tp\n")
    h.score(tmp_path, out)
    import json
    s = json.loads((out / "score.json").read_text())
    assert s["official_contexts"] == 30
    # TAZ (NaN, score 9) is dropped; C's 1e-14 is noise and scores 0 (amendment 1)
    assert s["mean_auroc_official"]["deltasignal"] == pytest.approx(1.0)
    assert s["mean_auroc_official"]["reachability"] == pytest.approx(0.75)
    assert s["mean_auroc_official"]["constant"] == pytest.approx(0.5)


def test_rank_min_matches_dreamtools_ties():
    # pandas rank(ascending=False, method="min"): ties share the best of their ranks.
    assert h.rank_min([0.7, 0.9, 0.7, 0.5]) == [2, 1, 2, 4]


def test_standings_and_place():
    table = {"a": {"c1": 0.9, "c2": 0.6}, "b": {"c1": 0.8, "c2": 0.8}, "c": {"c1": 0.5, "c2": 0.5}}
    auc, rk = h.standings(table, ["c1", "c2"])
    assert rk == {"a": 1.5, "b": 1.5, "c": 3.0}
    assert abs(auc["b"] - 0.8) < 1e-12
    assert h.place(auc, "a", True) == 2 and h.place(rk, "a", False) == 1
    # a resample repeating a context weights it twice
    _, rk2 = h.standings(table, ["c1", "c1"])
    assert rk2["a"] == 1.0


def test_rank_refuses_when_validation_fails(tmp_path):
    ctxs = [f"{c}_{s}" for c in h.CELL_LINES for s in h.OFFICIAL[c]]
    aucs = {ctx: {"t1": 0.6, "t2": 0.7, "t3": 0.5} for ctx in ctxs}
    (tmp_path / "aucs.json").write_text(json.dumps(aucs))
    rows = "Final rank,Team name,mean AUROC,mean Rank\n2,t1,0.9,2\n1,t2,0.9,1\n3,t3,0.9,3\n"
    (tmp_path / "lb.csv").write_text(rows)
    per = [{"cell": c, "stimulus": s, "deltasignal": 0.65, "reachability": 0.5}
           for c in h.CELL_LINES for s in h.OFFICIAL[c]]
    (tmp_path / "score.json").write_text(json.dumps({"per_context": per}))
    h.rank(tmp_path / "aucs.json", tmp_path / "lb.csv", tmp_path, boot=10)
    r = json.loads((tmp_path / "rank.json").read_text())
    assert len(r["validation_failures"]) == 3 and "deltasignal" not in r
    # with a matching sheet the entry is placed between t2 and t1
    # the sheet truncates AUROC and its ranks may run whole thirtieths high (amendment 3)
    (tmp_path / "lb.csv").write_text("Final rank,Team name,mean AUROC,mean Rank\n2,t1,0.6,2.1\n1,t2,0.7,1\n3,t3,0.5,3\n")
    h.rank(tmp_path / "aucs.json", tmp_path / "lb.csv", tmp_path, boot=10)
    r = json.loads((tmp_path / "rank.json").read_text())
    assert r["validation_failures"] == []
    assert r["deltasignal"]["place_by_mean_rank"] == 2 and r["deltasignal"]["place_by_mean_auroc_95"] == [2.0, 2.0]


def test_rank_validation_amendment_3(tmp_path):
    ctxs = [f"{c}_{s}" for c in h.CELL_LINES for s in h.OFFICIAL[c]]
    (tmp_path / "aucs.json").write_text(json.dumps({ctx: {"t1": 0.61239, "ML_BHK_KG": 0.7} for ctx in ctxs}))
    per = [{"cell": c, "stimulus": s, "deltasignal": 0.65, "reachability": 0.5}
           for c in h.CELL_LINES for s in h.OFFICIAL[c]]
    (tmp_path / "score.json").write_text(json.dumps({"per_context": per}))

    def fails(rows):
        (tmp_path / "lb.csv").write_text("Final rank,Team name,mean AUROC,mean Rank\n" + rows)
        h.rank(tmp_path / "aucs.json", tmp_path / "lb.csv", tmp_path, boot=2)
        return json.loads((tmp_path / "rank.json").read_text())["validation_failures"]
    assert fails("2,t1,0.6123,2\n1,MLBHKKG,0.7,1\n") == []            # truncated, alias
    assert len(fails("2,t1,0.6124,2\n1,MLBHKKG,0.7,1\n")) == 1         # rounded is not truncated
    assert len(fails("2,t1,0.6123,1.95\n1,MLBHKKG,0.7,1\n")) == 1      # sheet rank below ours
    assert len(fails("2,t1,0.6123,2.05\n1,MLBHKKG,0.7,1\n")) == 1      # not a whole thirtieth
    assert len(fails("2,t1,0.6123,2.1333333\n1,MLBHKKG,0.7,1\n")) == 1  # four thirtieths
    assert len(fails("1,t1,0.6123,2\n2,MLBHKKG,0.7,1\n")) == 1         # order differs
