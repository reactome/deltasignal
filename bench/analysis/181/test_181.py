import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ks_181 as K  # noqa: E402

SITES = {"p-X": {("X", 10)}, "p-Y": {("Y", 5)}}
KS = {"KIN": {("X", 10, "up"), ("Y", 7, "up"), ("Z", 3, "up")}, "PPASE": {("X", 10, "down")},
      "BOTH": {("X", 10, "up"), ("X", 10, "down")}, "OUTSIDE": {("X", 10, "up")}}
LUN = {"KIN", "PPASE", "BOTH", "Z"}


def test_readout_sites_parse_the_shipped_map():
    s = K.readout_sites()
    assert ("STAT3", 705) in s["p-STAT3"] and "beta-catenin" not in s


def test_arm_calls():
    calls = K.arm_calls(KS, SITES, {"Z": {"p-Y"}}, LUN)
    assert calls["K1"] == {("KIN", "p-X"): "up", ("PPASE", "p-X"): "down"}      # BOTH conflicts, OUTSIDE is not a Lun protein
    assert calls["K1p"][("KIN", "p-Y")] == "up"                                 # another site of the readout's protein
    assert calls["K2"] == {("KIN", "p-Y"): "up"}                                # KIN -> Z, Z reaches p-Y


def test_added_only_where_no_curated_route():
    cases = [{"gene": "KIN", "readout": "p-X", "reach": "up", "measured": "up"},
             {"gene": "PPASE", "readout": "p-X", "reach": "not", "measured": "down"}]
    add = K.added(cases, {("KIN", "p-X"): "up", ("PPASE", "p-X"): "down"})
    assert [c["gene"] for c in add] == ["PPASE"] and K.precision(add) == 1.0


def test_relabel_keeps_types_and_outsiders():
    out = K.relabel(KS, LUN, random.Random(1))
    assert out["OUTSIDE"] == KS["OUTSIDE"] and out["PPASE"] == KS["PPASE"]       # the only pure phosphatase maps to itself
    assert sorted(map(sorted, (out[e] for e in ("KIN", "BOTH")))) == sorted(map(sorted, (KS[e] for e in ("KIN", "BOTH"))))


def test_route_outcome_and_band():
    import routes_181 as R
    assert R.outcome({"measured": "up", "model": "not"}) == "silent"
    assert R.outcome({"measured": "down", "model": "up"}) == "wrong direction"
    assert R.outcome({"measured": "not", "model": "up"}) == "false alarm"
    assert [R.band(x) for x in (None, 0, 1, 5, 6, 20, 21)] == ["no path found", "0", "1-5", "1-5", "6-10", "11-20", "21+"]


def test_distances_count_reaction_steps():
    import routes_181 as R

    class Net:
        fwd = {"a": [("r1", "pos", "input")], "r1": [("b", "pos", "output")], "b": [("r2", "pos", "catalyst")],
               "r2": [("c", "pos", "output")]}
        kind = {"a": "entity", "b": "entity", "c": "entity", "r1": "reaction", "r2": "reaction"}
    assert R.distances(Net, ["a"])["c"] == (4, 2)
