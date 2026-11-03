import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import combine_132 as cm  # noqa: E402


def case(m, p, route, cluster="g"):
    return {"cluster": cluster, "measured": m, "predicted": p, "route": route, "pathway": "P"}


def test_balanced_accuracy_is_mean_recall_over_present_classes():
    assert cm.balanced_accuracy([("up", "up"), ("up", "not"), ("not", "not")]) == (0.5 + 1) / 2


def test_fallback_keeps_the_control_where_it_routes():
    control = {1: case("up", "not", True), 2: case("up", "not", False), 3: case("up", "up", True)}
    arm = {1: case("up", "up", True), 2: case("up", "up", True), 4: case("up", "up", True)}
    f = cm.fallback(control, arm)
    assert sorted(f) == [1, 2]                                    # only keys in both
    assert f[1]["predicted"] == "not" and f[2]["predicted"] == "up"


def test_ceiling_counts_control_errors_without_a_route_that_the_arm_reaches():
    control = {1: case("up", "not", False), 2: case("up", "not", True), 3: case("up", "up", False),
               4: case("down", "not", False)}
    arm = {1: case("up", "not", True), 2: case("up", "up", True), 3: case("up", "up", True),
           4: case("down", "not", False)}
    assert cm.ceiling(control, arm) == (3, 1)


def test_report_gate_pools_datasets_and_needs_both_point_estimates():
    control = {"experimental": {i: case("up", "not", False, f"e{i}") for i in range(4)},
               "korkut": {i: case("down", "down", True, f"k{i}") for i in range(4)}}
    arm = {"experimental": {i: case("up", "up", True, f"e{i}") for i in range(4)},
           "korkut": {i: case("down", "down", True, f"k{i}") for i in range(4)}}
    r = cm.report(control, arm, n=50)
    assert r["rule"]["experimental"]["fixed"] == 4 and r["rule"]["experimental"]["ceiling"] == 4
    assert r["rule"]["gate"]["pooled_gain"] > 0 and r["rule"]["gate"]["passes"]
    # under F the arm is read only where the control has no route: here all experimental cases
    assert r["rule+F"]["experimental"]["fixed"] == 4
    # a loss on one dataset fails the gate whatever the pooled interval says
    arm["korkut"] = {i: case("down", "not", True, f"k{i}") for i in range(4)}
    assert not cm.report(control, arm, n=50)["rule"]["gate"]["passes"]
