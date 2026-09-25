"""Pins for self_contained_inhibitor_pairs (the DS_SKIP_SELF_INH arm of
specs/012). Its first version inverted a multi-key index, so its answer
depended on PYTHONHASHSEED and it pooled activators across variant nodes."""
import json
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from benchmark_vs_mpbiopath import self_contained_inhibitor_pairs  # noqa: E402


def _write(d: Path):
    # Two variant nodes of reaction R-R1: vA takes ligand L1, vB takes L2.
    # Inhibitor complex I:L1 contains L1 and negatively regulates BOTH variants.
    # It is self-contained on vA only. Every node carries several member_leaves,
    # which is what made the old inverse index hash-order dependent.
    (d / "nodes.csv").write_text(
        "uuid,node_kind,diagram_entity_id,compartment,member_leaves,source_sets,chosen_members\n"
        "vA,reaction,R-R1,,,,\n"
        "vB,reaction,R-R1,,,,\n"
        "l1,simple_entity,R-L1,,R-L1|R-X|R-Y,,\n"
        "l2,simple_entity,R-L2,,R-L2|R-X|R-Z,,\n"
        "inh,simple_complex,R-IL1,,R-I|R-L1|R-X|R-Y,,\n")
    (d / "containment.csv").write_text(
        "stable_id,contains_stable_id,reactome_release\n"
        "R-IL1,R-IL1,97\nR-IL1,R-I,97\nR-IL1,R-L1,97\n")
    (d / "logic_network.csv").write_text(
        "source_id,target_id,pos_neg,and_or,edge_type\n"
        "l1,vA,pos,and,input\n"
        "l2,vB,pos,and,input\n"
        "inh,vA,neg,or,regulator\n"
        "inh,vB,neg,or,regulator\n")


def test_flags_only_the_variant_whose_own_activator_is_contained(tmp_path):
    _write(tmp_path)
    assert self_contained_inhibitor_pairs(tmp_path) == {("inh", "vA")}


def test_answer_does_not_depend_on_the_hash_seed(tmp_path):
    _write(tmp_path)
    code = ("import sys, json; sys.path.insert(0, %r); "
            "from benchmark_vs_mpbiopath import self_contained_inhibitor_pairs as f; "
            "print(json.dumps(sorted(f(__import__('pathlib').Path(%r)))))"
            % (str(Path(__file__).parent), str(tmp_path)))
    outs = set()
    for seed in ("1", "2", "3", "4"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        outs.add(subprocess.run([sys.executable, "-c", code], env=env,
                                capture_output=True, text=True, check=True).stdout)
    assert len(outs) == 1
    assert json.loads(outs.pop()) == [["inh", "vA"]]
