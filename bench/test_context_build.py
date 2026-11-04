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


def _typed(b, name, nodes, edges):
    """A pathway whose edges carry type and reaction id: (src, dst, edge_type, reaction stId)."""
    d = b / name
    d.mkdir(parents=True)
    (d / "nodes.csv").write_text(NODES + "".join(f"{u},{k},{st},,,,\n" for u, k, st in nodes))
    (d / "stid_to_uuid_mapping.csv").write_text("uuid,stable_id\n" + "".join(f"{u},{st}\n" for u, _, st in nodes))
    (d / "logic_network.csv").write_text(NET + "".join(f"{a},{t},pos,and,{et},1,{rid}\n" for a, t, et, rid in edges))
    (d / "cofactors.csv").write_text("stable_id,molecule,chebi_id,name,in_network,reactome_release\n")


def _curated_fixture(tmp_path):
    b = tmp_path / "base"
    b.mkdir()
    (b / "BUILD.json").write_text(json.dumps({"build_id": "base", "status": "complete"}))
    # P: E (root) is an input of reaction R1 -> Y; F (root) is also in P. Q: G -> S1 -> E, and S2 -> F.
    _typed(b, "R-HSA-1", [("pe", "simple_entity", "R-E"), ("pr", "reaction", "R-R1"), ("py", "simple_entity", "R-Y"),
                          ("pf", "simple_entity", "R-F"), ("pc", "simple_entity", "R-C")],
           [("pe", "pr", "input", "R-R1"), ("pr", "py", "output", "R-R1"), ("pf", "pr", "input", "R-R1"),
            ("pc", "pr", "catalyst", "")])
    _typed(b, "R-HSA-2", [("qg", "simple_entity", "R-G"), ("qs", "reaction", "R-S1"), ("qe", "simple_entity", "R-E"),
                          ("qs2", "reaction", "R-S2"), ("qf", "simple_entity", "R-F"), ("qc", "simple_entity", "R-C")],
           [("qg", "qs", "input", "R-S1"), ("qs", "qe", "output", "R-S1"), ("qs2", "qf", "output", "R-S2"),
            ("qs", "qc", "output", "R-S1")])
    (tmp_path / "links.tsv").write_text("pathway\tr\ts\tq\tentity\trole\tl2\n"
                                        "R-HSA-1\tR-R1\tR-S1\tx\tR-E\tinput\tFalse\n"
                                        # a catalyst edge carries no reaction id; it points at the reaction node
                                        "R-HSA-1\tR-R1\tR-S1\tx\tR-C\tcatalyst\tFalse\n")
    (tmp_path / "s_q.tsv").write_text("s\tq_diagram\tname\nR-S1\tR-HSA-2\tQ\n")
    (tmp_path / "centres.tsv").write_text("id\tpathway_name\nR-HSA-1\tP\n")
    return b


def test_curated_bridges_only_the_handed_off_entity(tmp_path):
    b = _curated_fixture(tmp_path)
    bridges, census = cb.plan_curated(b, tmp_path / "links.tsv", tmp_path / "s_q.tsv", {"R-HSA-1"})
    # E is curated as handed from S1 to R1: bridged. F is in both pathways but no curated handoff names it: not.
    assert bridges[("R-HSA-1", "R-HSA-2")]["in"] == [("qc", "pc"), ("qe", "pe")]
    assert [r["status"] for r in census] == ["bridged", "bridged"]
    # identity matching would bridge F too; the curated mode never does
    _, ident, _ = cb.plan(b, hub=5)
    assert ("qf", "pf") in ident[("R-HSA-1", "R-HSA-2")]["in"]


def test_curated_l2_drops_pairs_outside_encapsulated_pathways(tmp_path):
    b = _curated_fixture(tmp_path)
    bridges, census = cb.plan_curated(b, tmp_path / "links.tsv", tmp_path / "s_q.tsv", {"R-HSA-1"}, l2_only=True)
    assert bridges == {} and census == []


def test_curated_build_is_feed_in_and_an_empty_link_file_builds_the_control(tmp_path):
    b = _curated_fixture(tmp_path)
    cb.main([str(b), str(tmp_path / "c1"), "--curated", str(tmp_path / "links.tsv"), str(tmp_path / "s_q.tsv"),
             "--centres", str(tmp_path / "centres.tsv")])
    e = _edges(tmp_path / "c1" / "R-HSA-1")
    assert ("qe", "pe", "context_bridge") in e and ("qg", "qs", "input") in e
    assert not any(t == "context_bridge" and s.startswith("p") for s, _, t in e)
    assert not (tmp_path / "c1" / "R-HSA-2").exists()            # only centres are built
    (tmp_path / "empty.tsv").write_text("pathway\tr\ts\tq\tentity\trole\tl2\n")
    cb.main([str(b), str(tmp_path / "c0"), "--curated", str(tmp_path / "empty.tsv"), str(tmp_path / "s_q.tsv"),
             "--centres", str(tmp_path / "centres.tsv")])
    assert sorted(_edges(tmp_path / "c0" / "R-HSA-1")) == sorted(_edges(b / "R-HSA-1"))
    # every table the solver reads travels with the bundle (specs/100's RNA table was once left behind)
    (b / "R-HSA-1" / "rnas.csv").write_text("stable_id,schema_class,name,reactome_release\nR-E,EWAS,E mRNA,97\n")
    cb.main([str(b), str(tmp_path / "c0r"), "--curated", str(tmp_path / "empty.tsv"), str(tmp_path / "s_q.tsv"),
             "--centres", str(tmp_path / "centres.tsv")])
    assert (tmp_path / "c0r" / "R-HSA-1" / "rnas.csv").read_text() == (b / "R-HSA-1" / "rnas.csv").read_text()


def test_curated_set_handoff_bridges_each_member_to_the_same_member(tmp_path):
    b = _curated_fixture(tmp_path)
    # the handed-off entity R-SET is a set; the generator expanded it, so only its members have nodes
    with open(b / "R-HSA-1" / "nodes.csv", "a") as f:
        f.write("pm1,simple_entity,R-M1,,,,\n")
    with open(b / "R-HSA-1" / "logic_network.csv", "a") as f:
        f.write("pm1,pr,pos,and,input,1,R-R1\n")
    with open(b / "R-HSA-2" / "nodes.csv", "a") as f:
        f.write("qm1,simple_entity,R-M1,,,,\nqm2,simple_entity,R-M2,,,,\n")
    with open(b / "R-HSA-2" / "logic_network.csv", "a") as f:
        f.write("qs,qm1,pos,and,output,1,R-S1\nqs,qm2,pos,and,output,1,R-S1\n")
    (tmp_path / "set.tsv").write_text("pathway\tr\ts\tq\tentity\trole\tl2\nR-HSA-1\tR-R1\tR-S1\tx\tR-SET\tinput\tFalse\n")
    (tmp_path / "members.tsv").write_text("entity\tmember\tentity_class\nR-SET\tR-M1\tDefinedSet\nR-SET\tR-M2\tDefinedSet\n")
    bridges, census = cb.plan_curated(b, tmp_path / "set.tsv", tmp_path / "s_q.tsv", {"R-HSA-1"},
                                      members=tmp_path / "members.tsv")
    assert bridges[("R-HSA-1", "R-HSA-2")]["in"] == [("qm1", "pm1")]
    assert [r["status"] for r in census] == ["bridged"]
    # without the member table the set has no node in P
    _, census = cb.plan_curated(b, tmp_path / "set.tsv", tmp_path / "s_q.tsv", {"R-HSA-1"})
    assert [r["status"] for r in census] == ["entity_not_at_r_in_centre"]


def _transitive_fixture(tmp_path):
    b = _curated_fixture(tmp_path)
    # Q2 makes G (S0), which Q's S1 takes: a handoff into Q's root copy qg, two hops from P
    _typed(b, "R-HSA-3", [("ra", "simple_entity", "R-A"), ("rs0", "reaction", "R-S0"), ("rg", "simple_entity", "R-G")],
           [("ra", "rs0", "input", "R-S0"), ("rs0", "rg", "output", "R-S0")])
    head = "pathway\tr\ts\tq\tentity\trole\tl2\textracellular\n"
    (tmp_path / "tlinks.tsv").write_text(
        head + "R-HSA-1\tR-R1\tR-S1\tx\tR-E\tinput\tFalse\tFalse\n"
        "R-HSA-2\tR-S1\tR-S0\tx\tR-G\tinput\tFalse\tFalse\n"
        # P's R1 feeding back into Q: Q already joined above P's depth, so never bridged (no cycle)
        "R-HSA-2\tR-S1\tR-R1\tx\tR-G\tinput\tFalse\tFalse\n"
        # a cell-to-cell signal: listed, never bridged
        "R-HSA-1\tR-R1\tR-S2\tx\tR-F\tinput\tFalse\tTrue\n")
    (tmp_path / "ts_q.tsv").write_text("s\tq_diagram\tname\nR-S1\tR-HSA-2\tQ\nR-S0\tR-HSA-3\tQ2\nR-R1\tR-HSA-1\tP\n"
                                       "R-S2\tR-HSA-2\tQ\n")
    return b


def test_transitive_follows_handoffs_back_as_a_tree(tmp_path):
    b = _transitive_fixture(tmp_path)
    bridges, census = cb.plan_transitive(b, tmp_path / "tlinks.tsv", tmp_path / "ts_q.tsv", {"R-HSA-1"})
    assert bridges[("R-HSA-1", "R-HSA-2")]["in"] == [("qe", "pe")]
    assert bridges[("R-HSA-1", "R-HSA-3")]["in"] == [("rg", "qg")]          # depth 2, into Q's kept root
    status = {(r["pathway"], r["s"]): (r["status"], r["depth"]) for r in census}
    assert status[("R-HSA-1", "R-S1")] == ("bridged", 1)
    assert status[("R-HSA-2", "R-S0")] == ("bridged", 2)
    assert status[("R-HSA-2", "R-R1")][0] == "source_already_joined"
    assert status[("R-HSA-1", "R-S2")][0] == "extracellular"


def test_transitive_build_writes_every_hop_feed_in(tmp_path):
    b = _transitive_fixture(tmp_path)
    cb.main([str(b), str(tmp_path / "t"), "--curated", str(tmp_path / "tlinks.tsv"), str(tmp_path / "ts_q.tsv"),
             "--centres", str(tmp_path / "centres.tsv"), "--transitive"])
    e = _edges(tmp_path / "t" / "R-HSA-1")
    assert ("qe", "pe", "context_bridge") in e and ("rg", "qg", "context_bridge") in e
    assert ("ra", "rs0", "input") in e                                   # Q2's ancestors of its bridge source
    assert not any(t == "context_bridge" and s.startswith("p") for s, _, t in e)
    assert json.load(open(tmp_path / "t" / "BUILD.json"))["context_rule"] == "curated_transitive"
    # specs/133 amendment 1: every node's pathway, so the solver keeps identity rules inside one pathway
    origin = dict(line.split(",") for line in (tmp_path / "t" / "R-HSA-1" / "context_origin.csv").read_text().split()[1:])
    assert origin["pe"] == "R-HSA-1" and origin["qe"] == "R-HSA-2" and origin["rg"] == "R-HSA-3"
    assert {n for s, t, _ in e for n in (s, t)} <= set(origin)          # every node in the bundle has one


def test_whole_keeps_a_source_pathway_as_it_is_alone(tmp_path):
    # specs/133 amendment 2: Q also holds P's reaction R1 (its own copy qr) and a branch S2 -> F that never reaches P.
    # By default Q loses R1 (reduced) and S2 (pruned); under --whole it enters whole, so it solves as it does alone.
    b = _curated_fixture(tmp_path)
    q = b / "R-HSA-2"
    with open(q / "nodes.csv", "a") as f:
        f.write("qr,reaction,R-R1,,,,\nqy,simple_entity,R-Y,,,,\n")
    with open(q / "logic_network.csv", "a") as f:
        f.write("qe,qr,pos,and,input,1,R-R1\nqr,qy,pos,and,output,1,R-R1\n")
    for flag, out in (([], "d"), (["--whole"], "w")):
        cb.main([str(b), str(tmp_path / out), "--curated", str(tmp_path / "links.tsv"), str(tmp_path / "s_q.tsv"),
                 "--centres", str(tmp_path / "centres.tsv")] + flag)
    d, w = set(_edges(tmp_path / "d" / "R-HSA-1")), set(_edges(tmp_path / "w" / "R-HSA-1"))
    assert ("qe", "pe", "context_bridge") in d and ("qe", "pe", "context_bridge") in w      # same bridges
    assert ("qe", "qr", "input") not in d and ("qs2", "qf", "output") not in d
    assert {("qe", "qr", "input"), ("qr", "qy", "output"), ("qs2", "qf", "output")} <= w
    assert {e for e in w if e[2] != "context_bridge"} == set(_edges(b / "R-HSA-1")) | set(_edges(q))
    bj = json.load(open(tmp_path / "w" / "BUILD.json"))
    assert bj["context_whole"] is True and bj["context_pruned"] is False
