"""specs/167: the reached-base-rate reader."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis" / "167"))
import base_rate as b  # noqa: E402

unit = lambda k: k[0]  # noqa: E731

CASES = {("u1", "A"): ("down", "down", 0.5),   # reached, called, right
         ("u1", "B"): ("up", "down", 2.0),     # reached, called, a change but the wrong direction
         ("u1", "C"): ("down", "not", 0.2),    # reached, called, no change measured
         ("u1", "D"): ("not", "not", 0.9),     # reached, declined, no change
         ("u2", "E"): ("not", "up", 1.1),      # reached, declined, changed
         ("u2", "F"): ("not", "up", 1.0),      # unreached, changed
         ("u2", "G"): ("not", "not", 1.0),     # unreached
         ("u2", "H"): ("not", "not", 1.0)}     # unreached


def test_counts_and_rates():
    r = b.stats(CASES, unit, resamples=50)
    assert (r["cases"], r["units"], r["reached"], r["unreached"]) == (8, 2, 5, 3)
    assert r["base_rate_reached"] == 3 / 5 and r["base_rate_unreached"] == 1 / 3
    assert r["calls"] == 3 and r["call_rate"] == 3 / 5
    assert r["precision_change"] == 2 / 3 and r["precision_exact"] == 1 / 3 and r["sign_right"] == 1 / 2
    assert abs(r["lift"] - (2 / 3 - 3 / 5)) < 1e-12
    assert r["declined"] == 2 and r["declined_change_rate"] == 1 / 2


def test_calling_everything_reached_has_no_lift():
    cases = {("u", i): ("down", "down" if i < 3 else "not", 0.5) for i in range(10)}
    r = b.stats(cases, unit, resamples=50)
    assert r["call_rate"] == 1.0 and r["lift"] == 0.0 and r["lift_ci"] == [0.0, 0.0]
    assert r["declined"] == 0 and r["declined_change_rate"] is None


def test_an_informative_caller_has_positive_lift():
    cases = {}
    for u in range(20):
        for i in range(4):   # calls exactly the two readouts that changed, declines the two that did not
            cases[(u, i)] = ("down", "down", 0.5) if i < 2 else ("not", "not", 0.95)
    r = b.stats(cases, unit, resamples=200)
    assert r["precision_change"] == 1.0 and r["base_rate_reached"] == 0.5 and r["lift"] == 0.5
    assert r["lift_ci"][0] > 0 and r["declined_change_rate"] == 0.0


def test_bootstrap_resamples_units_not_cases():
    # one unit carries every call; resampling units makes the interval wide, resampling cases would not
    cases = {("big", i): ("down", "down", 0.5) for i in range(50)}
    for u in range(5):
        cases[(f"s{u}", 0)] = ("not", "not", 0.9)
    r = b.stats(cases, unit, resamples=500)
    assert r["units"] == 6 and r["lift_ci"][1] - r["lift_ci"][0] > 0.05


def test_no_reached_case_is_reported_as_absent():
    r = b.stats({("u", 0): ("not", "not", 1.0)}, unit, resamples=10)
    assert r["reached"] == 0 and r["lift"] is None and r["lift_ci"] == [None, None]
    assert r["base_rate_reached"] is None and r["precision_change"] is None


def test_decryptm_loader_keeps_scored_rows(tmp_path):
    p = tmp_path / "d.tsv"
    p.write_text("pathway_id\tdrug\tuniprot\tposition\tmeasured\tstatus\tpred_ui\tpredicted\n"
                 "P\tD\tQ1\t5\tdown\tscored\t0.5\tdown\n"
                 "P\tD\tQ1\t9\tnot\tno_readout\t\t\n")
    assert b.load_decryptm(p) == {("P", "D", "Q1", "5"): ("down", "down", 0.5)}
