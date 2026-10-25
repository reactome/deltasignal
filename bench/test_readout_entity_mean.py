"""specs/105: DS_KO_AGG=mean_entity weighs each curated entity a readout resolves to
once, not each variant copy."""
import os
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
from benchmark_vs_mpbiopath import entity_mean, load_output_entities  # noqa: E402


def test_each_entity_weighs_once():
    # ERBB2 readout: 1 EGFR copy knocked out beside 4 ERBB3 copies at baseline
    vals = {"e1": 0.0, "b1": 0.01, "b2": 0.01, "b3": 0.01, "b4": 0.01}
    groups = [["e1"], ["b1", "b2", "b3", "b4"]]
    assert entity_mean(vals, groups) == pytest.approx(0.005)
    # the plain mean counts copies: 0.008
    assert sum(vals.values()) / len(vals) == pytest.approx(0.008)


def test_one_group_is_the_plain_mean():
    vals = {"a": 0.02, "b": 0.01, "c": 0.0}
    assert entity_mean(vals, [["a", "b", "c"]]) == pytest.approx(0.01)
    assert entity_mean(vals, [["a", "b", "c"], []]) == pytest.approx(0.01)


def test_output_entities_skip_cofactors(tmp_path):
    (tmp_path / "nodes.csv").write_text(
        "uuid,node_kind,diagram_entity_id,compartment,member_leaves,source_sets,chosen_members\n"
        "r1,reaction,R-HSA-1963582,,,,\n"
        "o1,complex,R-HSA-1963593,,,,\n"
        "adp,small,R-ALL-29370,,,,\n")
    (tmp_path / "cofactors.csv").write_text(
        "stable_id,molecule,chebi_id,name,in_network,reactome_release\n"
        "R-ALL-29370,ADP,456216,ADP [cytosol],1,97\n")
    (tmp_path / "logic_network.csv").write_text(
        "source_id,target_id,pos_neg,and_or,edge_type,stoichiometry,edge_reaction_id\n"
        "r1,o1,pos,and,output,1,\n"
        "r1,adp,pos,and,output,1,\n"
        "x,r1,pos,and,input,1,\n")
    assert load_output_entities(tmp_path) == {"r1": {"R-HSA-1963593"}}


def test_unknown_mode_is_an_error():
    env = dict(os.environ, DS_KO_AGG="mean-entity")
    r = subprocess.run([sys.executable, "-c", "import benchmark_vs_mpbiopath"],
                       cwd=Path(__file__).parent, env=env, capture_output=True, text=True)
    assert r.returncode != 0 and "DS_KO_AGG must be one of" in r.stderr
