"""Pins for failure_structure.py on a network where every feature is known."""
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from failure_structure import Network, classify  # noqa: E402


def _net(tmp: Path, edges, nodes=None, containment=()):
    (tmp / "logic_network.csv").write_text(
        "source_id,target_id,pos_neg,and_or,edge_type\n"
        + "".join(f"{s},{t},{pn},and,{et}\n" for s, t, pn, et in edges))
    ids = sorted({x for s, t, _, _ in edges for x in (s, t)})
    (tmp / "nodes.csv").write_text(
        "uuid,node_kind,diagram_entity_id,compartment,member_leaves,source_sets,chosen_members\n"
        + "".join(f"{u},x,{(nodes or {}).get(u, 'R-' + u)},,,,\n" for u in ids))
    (tmp / "containment.csv").write_text(
        "stable_id,contains_stable_id,reactome_release\n"
        + "".join(f"{a},{b},97\n" for a, b in containment))
    return Network(tmp)


def case(pred, exp, d="0", g="A", o="Z"):
    return {"predicted": pred, "expected": exp, "direction": d, "gene_uuids": g, "output_uuids": o}


def test_sign_of_the_path_is_checked_against_the_expected_direction(tmp_path):
    net = _net(tmp_path, [("A", "B", "pos", "input"), ("B", "Z", "neg", "regulator")])
    # KD of A through one inhibition: Z should go UP. Expecting DOWN is wrong-sign only.
    assert classify(net, case("1", "2"))["reach"] == "signed_path"
    assert classify(net, case("1", "0"))["reach"] == "wrong_sign_only"
    assert classify(net, case("0", "1"))["error"] == "false_change"
    assert classify(net, case("1", "2"))["error"] == "missed"


def test_no_path(tmp_path):
    net = _net(tmp_path, [("A", "B", "pos", "input"), ("C", "Z", "pos", "input")])
    assert classify(net, case("1", "0"))["reach"] == "no_path"


def test_negative_and_positive_loops_and_welding(tmp_path):
    net = _net(tmp_path, [("A", "L1", "pos", "input"), ("L1", "L2", "pos", "input"),
                          ("L2", "L1", "neg", "regulator"), ("L2", "Z", "pos", "input")])
    r = classify(net, case("1", "0"))
    assert (r["loop"], r["loop_size"], r["derived_only"]) == ("has_inhibition", "small(<10)", "curated")
    net2 = _net(tmp_path, [("A", "L1", "pos", "input"), ("L1", "L2", "pos", "input"),
                           ("L2", "L1", "pos", "assembly"), ("L2", "Z", "pos", "input")])
    r2 = classify(net2, case("1", "0"))
    assert (r2["loop"], r2["derived_only"]) == ("all_positive", "welded")


def test_self_contained_inhibitor_on_the_route(tmp_path):
    net = _net(tmp_path, [("A", "R", "pos", "input"), ("A", "C", "pos", "assembly"),
                          ("C", "R", "neg", "regulator"), ("R", "Z", "pos", "output")],
               containment=[("R-C", "R-A")])
    r = classify(net, case("1", "0"))
    assert r["self_inh"] == "yes" and r["assembly"] == "yes"


from failure_structure import mh_risk_difference  # noqa: E402


def test_mh_removes_pure_between_pathway_composition():
    # P: all exposed, all wrong-prone; Q: all unexposed, all right -- but WITHIN
    # each pathway exposure makes no difference. Pooled differs; MH must be 0.
    rows = []
    for i in range(20):
        rows.append({"pathway": "P", "f": "y", "correct": i % 2 == 0})
        rows.append({"pathway": "P", "f": "n", "correct": i % 2 == 0})
        rows.append({"pathway": "Q", "f": "n", "correct": True})
    rd, used, n1, n0 = mh_risk_difference(rows, "f", {"y"}, {"n"})
    assert abs(rd) < 1e-12 and used == 1 and (n1, n0) == (20, 20)


def test_mh_recovers_a_real_within_pathway_effect():
    rows = [{"pathway": p, "f": "y", "correct": False} for p in "AB" for _ in range(10)]
    rows += [{"pathway": p, "f": "n", "correct": True} for p in "AB" for _ in range(10)]
    rd, used, _, _ = mh_risk_difference(rows, "f", {"y"}, {"n"})
    assert abs(rd - 1.0) < 1e-12 and used == 2


def test_self_inh_needs_a_path_from_the_shared_input(tmp_path):
    # C contains A but is not computed from A: the solver ignores it, so must this.
    net = _net(tmp_path, [("A", "R", "pos", "input"), ("B", "C", "pos", "assembly"),
                          ("C", "R", "neg", "regulator"), ("R", "Z", "pos", "output")],
               containment=[("R-C", "R-A")])
    assert net.self_inh == set()
