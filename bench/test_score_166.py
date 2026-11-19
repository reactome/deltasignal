"""specs/166: the literature-edge scorer."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis" / "166"))
import score_166 as s  # noqa: E402


def arm(calls):
    """calls: {key: (call, truth, fold)} on the korkut axis; experimental empty."""
    return {"korkut": dict(calls), "experimental": {}}


BASE = {("P", "d", "A"): ("not", "down", 1.0),     # unreached, wrong
        ("P", "d", "B"): ("down", "down", 0.5),    # reached, right
        ("P", "d", "C"): ("up", "not", 2.0),       # reached, wrong
        ("P", "d", "D"): ("not", "not", 1.0)}      # unreached, right


def test_three_numbers_split_new_from_previous():
    w = dict(BASE)
    w[("P", "d", "A")] = ("down", "down", 0.4)     # newly covered, fixed
    w[("P", "d", "D")] = ("up", "not", 1.5)        # newly covered, broken
    w[("P", "d", "C")] = ("not", "not", 1.1)       # previously covered, fixed
    keys = sorted(BASE)
    r = s.three_numbers(BASE, w, keys)
    assert r["coverage_gained"] == 2 and r["coverage_lost"] == 0
    assert r["new"]["accuracy"] == 0.5 and r["new"]["accuracy_without"] == 0.5
    assert r["new"]["change_calls"] == 2 and r["new"]["change_precision"] == 0.5
    assert r["previous"]["n"] == 2 and r["previous"]["fixed"] == 1 and r["previous"]["broken"] == 0
    assert r["all"]["net"] == 1


def test_coverage_lost_counts_a_reached_case_returning_to_one():
    w = dict(BASE)
    w[("P", "d", "B")] = ("not", "down", 1.0)
    r = s.three_numbers(BASE, w, sorted(BASE))
    assert r["coverage_lost"] == 1 and r["previous"]["broken"] == 1


def test_gate_needs_positive_net_and_every_seed_below():
    assert s.gate(3, [0, 1, 2] * 6 + [2])["pass"]           # 19 seeds, all below
    g = s.gate(3, [0, 3, 1])                       # a tie counts against
    assert not g["pass"] and g["seeds_at_or_above"] == 1 and g["rank_p"] == 0.5
    assert not s.gate(0, [-5, -4])["pass"]         # beats random, does not help
    assert s.gate(1, [0] * 19)["rank_p"] == 0.05


def test_gate_rank_p_is_reported_with_too_few_seeds():
    g = s.gate(5, [0, 0])
    assert g["beats_every_seed"] and g["rank_p"] == 1 / 3 and not g["pass"]


def test_score_uses_cases_common_to_every_primary_arm():
    w = dict(BASE)
    w[("P", "d", "A")] = ("down", "down", 0.4)
    sh = dict(BASE)
    del sh[("P", "d", "A")]                        # not scored in a shuffle: leaves the common set
    arms = {"without": arm(BASE), "with": arm(w), "shuf01": arm(sh)}
    res = s.score(arms)
    assert res["common"]["korkut"] == 3
    assert res["arms"]["with"]["primary_net"] == 0
    assert res["denominators"]["shuf01"]["korkut"] == 3


def test_ablation_arms_are_paired_with_without_only():
    w = dict(BASE)
    w[("P", "d", "A")] = ("down", "down", 0.4)
    sl = dict(w)
    arms = {"without": arm(BASE), "with": arm(w), "shuf01": arm(BASE), "only-x": arm(sl)}
    res = s.score(arms)
    assert res["arms"]["only-x"]["primary_net"] == 1
    assert res["gate"]["net_shuffles"] == [0]      # an ablation is not a seed


def test_unreached_gains_match_list_keys():
    w = dict(BASE)
    w[("P", "d", "A")] = ("down", "down", 0.4)
    keys = {"korkut": sorted(BASE), "curator_heldout": [], "curator_tuning": [], "experimental": []}
    base = {"korkut": BASE, "curator_heldout": {}, "curator_tuning": {}, "experimental": {}}
    a = dict(base, korkut=w)
    out = s.unreached_gains(base, a, keys, ({("P", "d", "A"), ("P", "d", "D")}, set()))
    assert out["korkut_list_scored"] == 2 and out["korkut_list_gained"] == 1
