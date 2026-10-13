"""specs/021: the decryptM drug bench's scoring, against hand-computed values."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import decryptm_bench as db  # noqa: E402


def test_readout_uuids_dedup_and_order():
    s2u = {"R-HSA-P1": ["a", "b"], "R-HSA-P2": ["b", "c"]}
    assert db.readout_uuids({"R-HSA-P2", "R-HSA-P1"}, s2u) == ["a", "b", "c"]
    assert db.readout_uuids(set(), s2u) == []


def test_score_case_excludes_pins_and_trivial_units():
    acts = {"r1": 0.002, "r2": 0.01, "p": 0.0}
    # pinned readout node left out of the mean: (0.2 + 1.0) / 2 = 0.6 -> DOWN
    st, ui, pred, rp = db.score_case(["r1", "r2", "p"], {"p"}, acts, {"p", "r1"})
    assert st == "scored" and abs(ui - 0.6) < 1e-12 and pred == db.bm.DOWN and rp == db.bm.DOWN
    # every readout node pinned: trivial, never scored as a DOWN hit
    assert db.score_case(["p"], {"p"}, acts, {"p"})[0] == "trivial"
    assert db.score_case([], {"p"}, acts, {"p"})[0] == "no_readout"
    assert db.score_case(["r1"], set(), acts, set())[0] == "no_pin"
    # unreached readout: reachability says NORMAL; a missing activity reads baseline
    st, ui, pred, rp = db.score_case(["zz"], {"p"}, acts, {"p"})
    assert (ui, pred, rp) == (1.0, db.bm.NORMAL, db.bm.NORMAL)


def test_mcnemar_exact():
    assert db.mcnemar_p(0, 0) == 1.0
    assert abs(db.mcnemar_p(5, 0) - 0.0625) < 1e-12        # 2 * (1/32)
    assert abs(db.mcnemar_p(1, 1) - 1.0) < 1e-12
    assert abs(db.mcnemar_p(8, 2) - 0.109375) < 1e-12      # 2 * 56/1024


def _row(drug, acc, pos, measured, pred, reach, status="scored", split="held-out"):
    return {"split": split, "status": status, "drug": drug, "uniprot": acc, "position": str(pos),
            "pathway": "P", "site_gene": "G", "measured": measured,
            "predicted": pred, "reach_predicted": reach}


def test_summarize_controls_and_baselines():
    rows = [
        _row("Dasatinib", "P27361", 202, "down", "down", "down"),     # ERK and dasatinib
        _row("SHP099", "P27361", 204, "down", "down", "down"),        # ERK only
        _row("Dasatinib", "Q1", 5, "up", "up", "down"),               # dasatinib only; DS-only hit
        _row("MK2206", "Q2", 7, "down", "not", "down"),               # neither; reach-only hit
        _row("MK2206", "Q3", 9, "not", "not", "down"),                # unchanged
        _row("MK2206", "Q4", 1, "down", "", "", status="trivial"),    # never scored
        _row("MK2206", "Q5", 1, "down", "down", "down", split="tuning"),
    ]
    s = db.summarize(rows)
    assert s["status"] == {"scored": 5, "trivial": 1}
    r = s["regulated"]
    assert (r["n"], r["deltasignal"], r["reachability"], r["always_down"]) == (4, 3, 3, 3)
    assert (r["ds_only"], r["reach_only"]) == (1, 1)
    assert s["regulated_no_erk"]["n"] == 2
    assert s["regulated_no_dasatinib"]["n"] == 2
    n = s["regulated_neither"]
    assert (n["n"], n["deltasignal"], n["reachability"]) == (1, 0, 1)
    assert (s["up"]["n"], s["up"]["deltasignal"], s["up"]["reachability"]) == (1, 1, 0)
    u = s["unchanged"]
    assert (u["n"], u["deltasignal"], u["reachability"], u["always_down"]) == (1, 1, 0, 0)
    assert set(s["by_drug"]) == {"Dasatinib", "SHP099", "MK2206"}
    # macro-F1 over the 5 scored cases. always-DOWN: F1(down) = 2*.6*1/1.6 = .75, others 0
    assert abs(s["macro_f1"]["always_down"] - 0.25) < 1e-12
