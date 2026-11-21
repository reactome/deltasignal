"""specs/168: the single-copy literature scorer and arm builder."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis" / "168"))
import score_168 as s  # noqa: E402
import build_arms_168 as b  # noqa: E402

BASE = {"korkut": {("P", "d", "A"): ("not", "down", 1.0),
                   ("P", "d", "B"): ("down", "down", 0.5),
                   ("P", "d", "C"): ("up", "not", 2.0)},
        "experimental": {("P", "g", "x", "r"): ("up", "up", 3.0)}}


def test_compare_counts_net_and_change_calls():
    arm = {"korkut": dict(BASE["korkut"]), "experimental": dict(BASE["experimental"])}
    arm["korkut"][("P", "d", "A")] = ("down", "down", 0.4)    # newly covered, fixed
    arm["korkut"][("P", "d", "B")] = ("not", "down", 0.9)     # previously covered, broken
    arm["korkut"][("P", "d", "C")] = ("not", "not", 1.1)      # previously covered, fixed
    r = s.compare(BASE, arm)
    assert r["korkut"]["all"] == {"fixed": 2, "broken": 1, "net": 1, "accuracy_without": 1 / 3,
                                  "accuracy_with": 2 / 3}
    assert r["korkut"]["change_calls_base"] == {"n": 2, "right": 1}
    assert r["korkut"]["change_calls_arm"] == {"n": 1, "right": 1}
    assert r["experimental"]["all"]["net"] == 0 and r["primary_net"] == 1


def test_compare_uses_only_cases_in_both_arms():
    arm = {"korkut": {("P", "d", "A"): ("down", "down", 0.4)}, "experimental": {}}
    r = s.compare(BASE, arm)
    assert r["korkut"]["cases"] == 1 and r["experimental"]["cases"] == 0


def test_gate_skips_shuffles_unless_net_positive():
    assert s.gate(0, [])["pass"] is False and "not run" in s.gate(-5, [])["shuffles"]
    g = s.gate(10, [1, 2, 3] * 7)
    assert g["pass"] and g["seeds_at_or_above"] == 0
    assert s.gate(10, [11] + [0] * 18)["pass"] is False


def test_builder_refuses_to_overwrite(tmp_path):
    src = tmp_path / "src"
    (src / "P").mkdir(parents=True)
    (tmp_path / "builds" / "X-with").mkdir(parents=True)
    with pytest.raises(SystemExit):
        b.build(src, tmp_path / "builds", "X", "with")
