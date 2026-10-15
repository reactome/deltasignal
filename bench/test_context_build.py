"""specs/056: the Reactome-context builder and how the bench reads its output."""
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import context_build as cb  # noqa: E402

NODES = "uuid,node_kind,diagram_entity_id,compartment,member_leaves,source_sets,chosen_members\n"
NET = "source_id,target_id,pos_neg,and_or,edge_type,stoichiometry,edge_reaction_id\n"


def _pathway(b, name, nodes, edges, cofactors=(), pools=""):
    d = b / name
    d.mkdir(parents=True)
    (d / "nodes.csv").write_text(NODES + "".join(f"{u},{k},{st},,,,\n" for u, k, st in nodes))
    (d / "stid_to_uuid_mapping.csv").write_text("uuid,stable_id\n" + "".join(f"{u},{st}\n" for u, _, st in nodes))
    (d / "logic_network.csv").write_text(NET + "".join(f"{a},{t},pos,and,input,1,\n" for a, t in edges))
    (d / "cofactors.csv").write_text("stable_id,molecule,chebi_id,name,in_network,reactome_release\n"
                                     + "".join(f"{c},x,1,x,1,97\n" for c in cofactors))
    (d / "pools.csv").write_text("pool_id,node_uuid,stable_id,role,is_base\n" + pools)


def _edges(d):
    with open(d / "logic_network.csv", newline="") as f:
        return [(r["source_id"], r["target_id"], r["edge_type"]) for r in csv.DictReader(f)]


def _fixture(tmp_path):
    b = tmp_path / "base"
    b.mkdir()
    (b / "BUILD.json").write_text(json.dumps({"build_id": "base", "status": "complete"}))
    # P: X is a root that feeds readout Y. Q: U -> r -> X (X produced in Q), plus
    # an unrelated branch Z -> W that cannot reach P, and ATP produced in Q.
    _pathway(b, "R-HSA-1", [("px", "simple_entity", "R-X"), ("py", "simple_entity", "R-Y"),
                            ("patp", "simple_entity", "R-ATP")],
             [("px", "py"), ("patp", "py")], cofactors=["R-ATP"])
    _pathway(b, "R-HSA-2", [("qu", "simple_entity", "R-U"), ("qr", "reaction", "R-R"),
                            ("qx", "simple_entity", "R-X"), ("qz", "simple_entity", "R-Z"),
                            ("qw", "simple_entity", "R-W"), ("qatp", "simple_entity", "R-ATP"),
                            ("qa", "simple_entity", "R-ATPIN")],
             [("qu", "qr"), ("qr", "qx"), ("qz", "qw"), ("qa", "qatp")], cofactors=["R-ATP"],
             pools="pool1,qz,R-Z,state,True\npool1,qw,R-W,state,False\n")
    return b


def test_bridge_runs_from_producer_to_root_and_skips_cofactors(tmp_path):
    b = _fixture(tmp_path)
    names, bridges, _ = cb.plan(b, hub=5)
    assert bridges[("R-HSA-1", "R-HSA-2")]["in"] == [("qx", "px")]
    assert bridges[("R-HSA-1", "R-HSA-2")]["out"] == []    # px is not produced in P
    # seen from Q, P only receives (px is a root in P), so P is no neighbour of Q
    assert bridges[("R-HSA-2", "R-HSA-1")]["in"] == []
    assert all("patp" not in (s, t) for bs in bridges.values() for s, t in bs["in"] + bs["out"])


def test_hub_entities_are_never_bridged(tmp_path):
    b = _fixture(tmp_path)
    _, bridges, n_paths = cb.plan(b, hub=1)          # R-X is in 2 pathways: a hub at cutoff 1
    assert n_paths["R-X"] == 2
    assert bridges == {}


def test_context_bundle_prunes_and_scores_own_nodes(tmp_path):
    b = _fixture(tmp_path)
    out = tmp_path / "ctx"
    cb.main([str(b), str(out)])
    d = out / "R-HSA-1"
    e = _edges(d)
    assert ("qx", "px", "context_bridge") in e
    assert ("qu", "qr", "input") in e and ("qr", "qx", "input") in e
    # Z -> W cannot reach P, and its pool goes with it
    assert not any("qz" in (s, t) for s, t, _ in e)
    with open(d / "pools.csv", newline="") as f:
        assert list(csv.DictReader(f)) == []
    with open(d / "context_own.csv", newline="") as f:
        assert {r["uuid"] for r in csv.DictReader(f)} == {"px", "py", "patp"}
    with open(d / "nodes.csv", newline="") as f:
        assert {r["uuid"] for r in csv.DictReader(f)} == {"px", "py", "patp", "qu", "qr", "qx"}
    meta = json.loads((out / "BUILD.json").read_text())
    assert meta["context_of"] == "base" and meta["context_hub"] == 5 and meta["status"] == "complete"
    # Q's own bundle: nothing P produces is a root in Q, so nothing enters Q
    assert meta["context"]["R-HSA-2"]["neighbours"] == []
    assert not any(t == "context_bridge" for _, _, t in _edges(out / "R-HSA-2"))


def test_pool_ids_are_prefixed_by_pathway(tmp_path):
    b = _fixture(tmp_path)
    out = tmp_path / "ctx"
    cb.main([str(b), str(out), "--no-prune"])
    with open(out / "R-HSA-1" / "pools.csv", newline="") as f:
        assert {r["pool_id"] for r in csv.DictReader(f)} == {"R-HSA-2:pool1"}


def test_bench_reads_roots_without_bridges_and_scores_own_readouts(tmp_path, monkeypatch):
    import benchmark_vs_mpbiopath as bm
    b = _fixture(tmp_path)
    out = tmp_path / "ctx"
    cb.main([str(b), str(out)])
    monkeypatch.setattr(bm, "PIN_SCOPE", "root")
    monkeypatch.setattr(bm, "PIN_NETWORK_ROOT", "")
    # the gene's root forms: px (a root in P, bridged from Q) and qu (a root in Q)
    pins, indeg = bm.pin_occurrences(out / "R-HSA-1", {"G": ["px", "qu", "qx"]}, {})
    assert pins == {"G": ["px", "qu"]}
    assert bm.context_own_uuids(out / "R-HSA-1") == {"px", "py", "patp"}
    assert bm.context_own_uuids(b / "R-HSA-1") is None
    monkeypatch.setattr(bm, "CATALOG_ROOT", out)
    assert bm.context_label() == "base+hub5"
    (out / "BUILD.json").write_text(json.dumps({"context_of": "base", "context_hub": 5, "context_feedback": False}))
    assert bm.context_label() == "base+hub5,feed-in"
    monkeypatch.setattr(bm, "CATALOG_ROOT", b)
    assert bm.context_label() == "-"


def test_a_neighbour_loses_the_centres_reactions(tmp_path):
    # P has reaction R-R (u -> r -> x) and Q is its parent: the same reaction
    # plus one P lacks (x -> s -> y). Q's copy of R-R is dropped, so Q's x is a
    # root there and is fed from P's produced x; P's own copy of R-R is never
    # re-linked to Q's.
    b = tmp_path / "base"
    b.mkdir()
    (b / "BUILD.json").write_text(json.dumps({"build_id": "base", "status": "complete"}))
    _pathway(b, "R-HSA-1", [("pu", "simple_entity", "R-U"), ("pr", "reaction", "R-R"), ("px", "simple_entity", "R-X"),
                            ("py", "simple_entity", "R-Y")], [("pu", "pr"), ("pr", "px"), ("py", "pu")])
    _pathway(b, "R-HSA-2", [("qu", "simple_entity", "R-U"), ("qr", "reaction", "R-R"), ("qx", "simple_entity", "R-X"),
                            ("qs", "reaction", "R-S"), ("qy", "simple_entity", "R-Y")],
             [("qu", "qr"), ("qr", "qx"), ("qx", "qs"), ("qs", "qy")])
    _, bridges, _ = cb.plan(b, hub=5)
    bq = bridges[("R-HSA-1", "R-HSA-2")]
    assert bq["out"] == [("pu", "qu"), ("px", "qx")] and bq["in"] == [("qy", "py")]
    assert all("qr" not in (a, c) for a, c, _ in bq["edges"])
    out = tmp_path / "ctx"
    cb.main([str(b), str(out)])
    e = _edges(out / "R-HSA-1")
    assert ("px", "qx", "context_bridge") in e and ("qy", "py", "context_bridge") in e
    assert not any("qr" in (a, c) or "qu" in (a, c) for a, c, _ in e)   # qu fed only Q's copy of R-R


def test_feed_in_adds_no_bridge_out_of_the_centre(tmp_path):
    b = tmp_path / "base"
    b.mkdir()
    (b / "BUILD.json").write_text(json.dumps({"build_id": "base", "status": "complete"}))
    # P: u -> x, and y a root. Q: x -> s -> y. Feedback: P's x -> Q's x -> y -> P's y.
    _pathway(b, "R-HSA-1", [("pu", "simple_entity", "R-U"), ("px", "simple_entity", "R-X"),
                            ("py", "simple_entity", "R-Y")], [("pu", "px"), ("py", "pu")])
    _pathway(b, "R-HSA-2", [("qx", "simple_entity", "R-X"), ("qs", "reaction", "R-S"),
                            ("qy", "simple_entity", "R-Y")], [("qx", "qs"), ("qs", "qy")])
    cb.main([str(b), str(tmp_path / "fb")])
    cb.main([str(b), str(tmp_path / "in"), "--feed-in"])
    fb, fi = _edges(tmp_path / "fb" / "R-HSA-1"), _edges(tmp_path / "in" / "R-HSA-1")
    assert ("px", "qx", "context_bridge") in fb and ("qy", "py", "context_bridge") in fb
    assert ("qy", "py", "context_bridge") in fi and ("px", "qx", "context_bridge") not in fi
    assert json.loads((tmp_path / "in" / "BUILD.json").read_text())["context_feedback"] is False
