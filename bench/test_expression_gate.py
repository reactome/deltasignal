import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import expression_gate as g  # noqa: E402

BELOW = {"SET1": {"A", "B"}, "CPLX": {"B", "C"}}


def test_complex_needs_every_leaf():
    assert g.required_leaves({"node_kind": "simple_complex", "member_leaves": "X|Y"}, BELOW) == ("any", {"X", "Y"})


def test_set_variant_needs_fixed_leaves_and_its_chosen_member_only():
    row = {"node_kind": "set_variant", "member_leaves": "X|A|B", "source_sets": "SET1", "chosen_members": "A"}
    assert g.required_leaves(row, BELOW) == ("any", {"X", "A"})            # B is the member not chosen


def test_set_variant_with_a_complex_member_needs_its_leaves():
    row = {"node_kind": "set_variant", "member_leaves": "X|B|C|A", "source_sets": "SET1", "chosen_members": "CPLX"}
    assert g.required_leaves(row, {**BELOW, "SET1": {"A", "B", "C"}}) == ("any", {"X", "B", "C"})


def test_pools_and_plain_sets():
    assert g.required_leaves({"node_kind": "set_pool", "member_leaves": "A|B"}, BELOW) == ("any", set())
    assert g.required_leaves({"node_kind": "simple_entity", "member_leaves": "A|B"}, BELOW) == ("all", {"A", "B"})


def test_gated_nodes(tmp_path):
    (tmp_path / "nodes.csv").write_text(
        "uuid,node_kind,diagram_entity_id,compartment,member_leaves,source_sets,chosen_members\n"
        "c,simple_complex,,,X|Y,,\nv1,set_variant,,,X|A|B,SET1,A\nv2,set_variant,,,Z|A|B,SET1,B\n"
        "s,simple_entity,,,A|B,,\np,set_pool,,,A|B,,\n")
    assert g.gated_nodes(tmp_path, {"A"}, BELOW) == {"v1"}
    assert g.gated_nodes(tmp_path, {"A", "B"}, BELOW) == {"v1", "v2", "s"}


def test_absent_genes(tmp_path):
    f = tmp_path / "e.csv"
    f.write_text(",AKT1 (207),EGFR (1956),MYC (4609)\nACH-1,5.0,0.2,\nACH-2,0.1,3,1\n")
    assert g.absent_genes(f, "ACH-1", {"AKT1", "EGFR", "MYC", "NOPE"}) == {"EGFR"}
    assert g.absent_genes(f, "ACH-9", {"AKT1"}) is None
