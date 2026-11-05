"""specs/139: the target-activity gate's pieces (values-free fixtures)."""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis" / "139"))
import markers as mk  # noqa: E402
import score_139 as sc  # noqa: E402

UP = math.log2(2)
DOWN = -math.log2(2)


def test_ranks_within_set_mid_rank_and_minimum_lines():
    base = {("s1", "a"): {"m": 1.0}, ("s1", "b"): {"m": 2.0}, ("s1", "c"): {"m": 2.0}, ("s1", "d"): {"m": 3.0},
            ("s2", "a"): {"m": 9.0}, ("s2", "b"): {"m": 0.0}}
    r = sc.ranks(base, min_lines=4)
    assert r[("s1", "a")]["m"] == 0.0
    assert r[("s1", "b")]["m"] == r[("s1", "c")]["m"] == 0.5     # tie: mid-rank
    assert r[("s1", "d")]["m"] == 1.0
    assert ("s2", "a") not in r                                   # two lines: no rank


def test_activity_mean_over_measured_markers():
    table = {("s1", "a"): {"m1": 0.0, "m2": 1.0}}
    s = {"set": "s1", "cell_line": "a"}
    assert sc.activity(s, {"m1", "m2", "absent"}, table) == 0.5
    assert sc.activity(s, {"absent"}, table) is None
    assert sc.activity({"set": "s9", "cell_line": "a"}, {"m1"}, table) is None


def test_gate_turns_every_call_to_no_change():
    pairs = [(DOWN, DOWN), (DOWN, 0.0), (UP, UP), (0.0, 0.0)]
    tot, hit, gated = sc.tally(pairs)
    assert tot == {-1: 1, 0: 2, 1: 1}
    assert hit == {-1: 1, 0: 1, 1: 1}
    assert gated == {0: 2}                        # gated: only the measured no-changes are right
    assert sc.gated_at({}, 0.1, 0.2) and not sc.gated_at({}, 0.3, 0.2)
    assert not sc.gated_at({}, None, 0.6) and not sc.gated_at({}, 0.0, None)


def test_choose_prefers_identity_on_a_tie():
    recs = sc.build([{"drug": "D", "cell_line": "a", "values": {"x": 0.0}}], {("D", "x"): 0.0}, {"x"}, {})
    assert sc.choose(recs, [0.1], [0]) is None    # gating a sample already called NOT changes nothing


def test_build_leaves_markers_out():
    s = {"drug": "D", "cell_line": "a", "values": {"m": DOWN, "x": DOWN, "y": 0.0}}
    recs = sc.build([s], {("D", "m"): DOWN, ("D", "x"): DOWN, ("D", "y"): 0.0}, {"m", "x", "y"}, {"D": {"m"}})
    assert [a for a, _ in recs[0]["pairs"]] == [DOWN, 0.0]


def test_primary_targets_kinobeads_within_ten_fold():
    rows = [{"source": "kinobeads", "kd_app_nM": "10", "target": "P1"},
            {"source": "kinobeads", "kd_app_nM": "100", "target": "P2"},
            {"source": "kinobeads", "kd_app_nM": "101", "target": "P3"}]
    assert mk.primary_targets(rows, {}) == {"P1", "P2"}
    assert mk.primary_targets([{"source": "binding", "kd_app_nM": "", "target": "P4"}], {}) == set()
    assert mk.primary_targets([{"source": "class", "kd_app_nM": "", "target": "MAP2K1"}],
                              {"MAP2K1": {"Q02750"}}) == {"Q02750"}


def test_markers_need_a_target_catalyst():
    ab = {"ERK_p": ["MAPK3:202"], "MEK_p": ["MAP2K1:217"]}
    cat = {"MAPK3:202": {"Q02750"}, "MAP2K1:217": {"P15056"}}
    assert mk.markers(ab, cat, {"Q02750"}) == ["ERK_p"]
