import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import drug_activity as da  # noqa: E402

NODES = ("uuid,node_kind,diagram_entity_id,compartment,member_leaves,source_sets,chosen_members\n"
         "mek,simple_entity,MEK,,LMEK,,\nrafmek,simple_complex,RAFMEK,,LRAF|LMEK,,\nmekerk,simple_complex,MEKERK,,LMEK|LERK,,\n"
         "r1,reaction,R1,,,,\nr2,reaction,R2,,,,\nr3,reaction,R3,,,,\nerk,simple_entity,ERK,,LERK,,\n")
EDGES = ("source_id,target_id,pos_neg,and_or,edge_type,stoichiometry,edge_reaction_id\n"
         "mek,r1,pos,and,catalyst,1,R1\nrafmek,r2,pos,and,catalyst,1,R2\nmekerk,r3,pos,and,catalyst,1,R3\n"
         "mek,r3,pos,and,input,1,R3\nerk,r1,pos,and,input,1,R1\n")


def _net(tmp_path):
    (tmp_path / "nodes.csv").write_text(NODES)
    (tmp_path / "logic_network.csv").write_text(EDGES)


def test_cut_respects_active_units(tmp_path):
    _net(tmp_path)
    units = {("R1", "MEK"): None,                      # no active unit: the holder's activity is cut
             ("R2", "RAFMEK"): {"LRAF"},               # RAF is the active unit: a MEK drug leaves it
             ("R3", "MEKERK"): {"LMEK"}}               # MEK is the active unit
    assert da.cut_edges(tmp_path, {"LMEK"}, units) == {("mek", "r1"), ("mekerk", "r3")}   # inputs never cut


def test_redirect_keeps_the_edge_and_holds_its_source_at_zero():
    parsed = {"nodes": [{"uuid": "mek", "baseline": 0.01, "entity_type": "x"}, {"uuid": "r1", "baseline": 0.01}],
              "edges": [{"parent_uuid": "mek", "child_uuid": "r1", "edge_type": "catalyst"},
                        {"parent_uuid": "mek", "child_uuid": "r9", "edge_type": "input"}],
              "pathways": []}
    payload, obs = da.redirected(parsed, {("mek", "r1")}, {})
    assert [e["parent_uuid"] for e in payload["edges"]] == ["drugblock:mek:r1", "mek"]
    assert obs == {"drugblock:mek:r1": [0.0, 1.0]}
    assert len(payload["nodes"]) == 3 and parsed["edges"][0]["parent_uuid"] == "mek"     # the original is untouched


TABLE = "# test\naxis\tdrug\tlabel\nkorkut\tAK\tknockout\nkorkut\t901\tactivity\nkorkut\tRY\tunknown\ncppa\tAK\tactivity\n"


def test_mechanism_labels_per_axis(tmp_path, monkeypatch):
    (tmp_path / "t.tsv").write_text(TABLE)
    monkeypatch.setattr(da, "DRUG_ACTION", "mechanism")
    monkeypatch.setattr(da, "MECHANISM_UNKNOWN", "knockout")
    assert da.drug_actions("korkut", ["AK", "901", "RY"], tmp_path / "t.tsv") == \
        {"AK": "knockout", "901": "activity", "RY": "knockout"}
    monkeypatch.setattr(da, "MECHANISM_UNKNOWN", "activity")
    assert da.drug_actions("korkut", ["RY"], tmp_path / "t.tsv") == {"RY": "activity"}
    assert da.drug_actions("cppa", ["AK"], tmp_path / "t.tsv") == {"AK": "activity"}     # rows of its own axis


def test_mechanism_refuses_an_unlabelled_drug(tmp_path, monkeypatch):
    (tmp_path / "t.tsv").write_text(TABLE)
    monkeypatch.setattr(da, "DRUG_ACTION", "mechanism")
    import pytest
    with pytest.raises(SystemExit, match="no korkut label"):
        da.drug_actions("korkut", ["AK", "PLX"], tmp_path / "t.tsv")


def test_other_modes_ignore_the_table(monkeypatch):
    monkeypatch.setattr(da, "DRUG_ACTION", "activity")
    assert da.drug_actions("korkut", ["AK"], "/nonexistent") == {"AK": "activity"}


def test_knockout_wins_a_gene_two_members_share():
    assert da.gene_actions([("activity", ["MAP2K1", "AKT1"]), ("knockout", ["AKT1"])]) == \
        {"MAP2K1": "activity", "AKT1": "knockout"}


def test_mixed_perturbation_pins_knockout_genes_and_cuts_activity_genes(tmp_path):
    _net(tmp_path)

    class BM:
        PERTURB_UI_DOWN, PIN_CONFIDENCE = 0.0, 1.0

        @staticmethod
        def pin_occurrences(pdir, g2u, gene_to_stids):
            return {g: us for g, us in g2u.items() if us}, None

    g2s = {"MEK1": ["LMEK"], "ERK1": ["LERK"]}
    s2u = {"LMEK": ["mek"], "LERK": ["erk"]}
    units = {("R1", "MEK"): None, ("R2", "RAFMEK"): {"LRAF"}, ("R3", "MEKERK"): {"LMEK"}}
    pins, obs, cut, roots = da.perturb(BM, tmp_path, None, ["ERK1", "MEK1"], g2s, s2u, units,
                                       {"MEK1": "activity", "ERK1": "knockout"})
    assert obs == {"erk": [0.0, 1.0]} and cut == {("mek", "r1"), ("mekerk", "r3")} and roots == {"mek", "erk"}
    # a single action string reproduces specs/122's two modes
    assert da.perturb(BM, tmp_path, None, ["MEK1"], g2s, s2u, units, "knockout")[1:3] == ({"mek": [0.0, 1.0]}, set())
    assert da.perturb(BM, tmp_path, None, ["MEK1"], g2s, s2u, units, "activity")[1] == {}


def test_the_committed_table_labels_every_korkut_drug_and_hpn_dream(monkeypatch):
    monkeypatch.setattr(da, "DRUG_ACTION", "mechanism")
    import csv as _csv
    codes = [r["code"] for r in _csv.DictReader((l for l in open(Path(__file__).parent / "analysis/062/drugs.tsv")
                                                  if not l.startswith("#")), delimiter="\t")]
    assert set(da.drug_actions("korkut", codes).values()) <= {"knockout", "activity"}
    assert da.drug_actions("hpn_dream", ["AZD8055"]) == {"AZD8055": "activity"}
