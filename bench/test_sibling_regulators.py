"""Pins for sibling_regulator_pairs (specs/025): a regulator naming a specific
set member regulates only the variants built from that member."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from benchmark_vs_mpbiopath import sibling_regulator_pairs  # noqa: E402


def _write(d: Path):
    # Reaction R-R1 has variants va (CREBBP + IRF3) and vb (EP300 + IRF3).
    # CREBBP:NS1 negatively regulates both; it names CREBBP.
    (d / "nodes.csv").write_text(
        "uuid,node_kind,diagram_entity_id,compartment,member_leaves,source_sets,chosen_members\n"
        "va,reaction,R-R1,,,,\n"
        "vb,reaction,R-R1,,,,\n"
        "crebbp,simple_entity,R-CREBBP,,,,\n"
        "ep300,simple_entity,R-EP300,,,,\n"
        "irf3,simple_entity,R-IRF3,,,,\n"
        "ns1c,simple_complex,R-CNS1,,R-CREBBP|R-NS1,,\n"
        "other,simple_complex,R-OTHER,,R-X|R-Y,,\n")
    (d / "logic_network.csv").write_text(
        "source_id,target_id,pos_neg,and_or,edge_type\n"
        "crebbp,va,pos,and,input\nirf3,va,pos,and,input\n"
        "ep300,vb,pos,and,input\nirf3,vb,pos,and,input\n"
        "ns1c,va,neg,or,regulator\nns1c,vb,neg,or,regulator\n"
        "other,va,neg,or,regulator\nother,vb,neg,or,regulator\n")


def test_member_specific_regulator_leaves_only_the_sibling_variant(tmp_path):
    _write(tmp_path)
    assert sibling_regulator_pairs(tmp_path) == {("ns1c", "vb")}


def test_a_regulator_naming_no_set_member_stays_on_every_variant(tmp_path):
    _write(tmp_path)
    assert ("other", "va") not in sibling_regulator_pairs(tmp_path)
    assert ("other", "vb") not in sibling_regulator_pairs(tmp_path)


def test_a_reaction_without_variants_is_untouched(tmp_path):
    (tmp_path / "nodes.csv").write_text(
        "uuid,node_kind,diagram_entity_id,compartment,member_leaves,source_sets,chosen_members\n"
        "r,reaction,R-R,,,,\na,simple_entity,R-A,,,,\nc,simple_complex,R-C,,R-A|R-B,,\n")
    (tmp_path / "logic_network.csv").write_text(
        "source_id,target_id,pos_neg,and_or,edge_type\na,r,pos,and,input\nc,r,neg,or,regulator\n")
    assert sibling_regulator_pairs(tmp_path) == set()


def test_a_variant_with_no_inputs_is_left_alone(tmp_path):
    _write(tmp_path)
    # add a third variant of R-R1 with no inputs; its regulators must stay
    (tmp_path / "nodes.csv").write_text((tmp_path / "nodes.csv").read_text() + "vc,reaction,R-R1,,,,\n")
    (tmp_path / "logic_network.csv").write_text(
        (tmp_path / "logic_network.csv").read_text() + "ns1c,vc,neg,or,regulator\n")
    assert ("ns1c", "vc") not in sibling_regulator_pairs(tmp_path)
