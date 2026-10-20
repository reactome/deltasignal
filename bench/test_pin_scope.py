"""Pins for entry_occurrences (specs/023): pin where the gene enters the
network, compute everything downstream."""
import sys
import pytest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from benchmark_vs_mpbiopath import entry_occurrences  # noqa: E402


def test_root_entity_is_pinned_and_its_complexes_are_not():
    # A (the gene) + B -> AB complex -> reaction R -> AB:C
    adj = {"A": ["AB"], "B": ["AB"], "AB": ["R"], "R": ["ABC"]}
    assert entry_occurrences(["A", "AB", "ABC"], adj) == ["A"]


def test_independent_occurrences_are_each_an_entry():
    # Two branches that never meet: both are where the gene enters.
    adj = {"A1": ["X"], "A2": ["Y"]}
    assert entry_occurrences(["A1", "A2"], adj) == ["A1", "A2"]


def test_first_mid_pathway_occurrence_is_pinned_when_there_is_no_root():
    # The gene is only produced by a reaction (transcription) and flows on.
    adj = {"DNA": ["TX"], "TX": ["P"], "P": ["PC"]}
    assert entry_occurrences(["P", "PC"], adj) == ["P"]


def test_members_that_only_reach_each_other_fall_back_to_all():
    adj = {"A": ["B"], "B": ["A"]}
    assert sorted(entry_occurrences(["A", "B"], adj)) == ["A", "B"]


def test_a_cycle_downstream_of_the_root_is_not_pinned():
    adj = {"A": ["C1"], "C1": ["C2"], "C2": ["C1"]}
    assert entry_occurrences(["A", "C1", "C2"], adj) == ["A"]


from benchmark_vs_mpbiopath import root_occurrences  # noqa: E402


def test_root_scope_keeps_only_nodes_with_no_incoming_edge():
    indeg = {"AB": 2, "R": 1, "ABC": 1}           # A and B are roots
    assert root_occurrences(["A", "AB", "ABC"], indeg) == ["A"]


def test_a_root_complex_containing_the_gene_is_pinned():
    # A root complex that is not decomposed still counts: it contains the gene.
    assert root_occurrences(["AB", "ABC"], {"ABC": 1}) == ["AB"]


def test_root_scope_does_not_invent_a_pin_for_a_gene_with_no_root_form():
    # Unlike `entry`, a gene only produced mid-pathway is NOT perturbed.
    assert root_occurrences(["P", "PC"], {"P": 1, "PC": 1}) == []


def test_root_scope_drops_duplicates():
    assert root_occurrences(["A", "A", "B"], {}) == ["A", "B"]


from benchmark_vs_mpbiopath import recycled_root_occurrences  # noqa: E402


def _graph(edges):
    inc, fwd = {}, {}
    for s, t, et in edges:
        inc.setdefault(t, []).append((s, et))
        fwd.setdefault(s, []).append(t)
    return inc, fwd


def test_enzyme_regenerated_by_its_own_cycle_is_its_root():
    # E + S -> ES -> E + P : the free enzyme is fed only by its own cycle.
    inc, fwd = _graph([("E", "bind", "input"), ("S", "bind", "input"), ("bind", "ES", "output"),
                       ("ES", "cat", "input"), ("cat", "E", "output"), ("cat", "P", "output")])
    assert recycled_root_occurrences(["E", "ES"], {"E"}, inc, fwd) == ["E"]


def test_an_entity_produced_from_outside_is_not_a_recycled_root():
    inc, fwd = _graph([("X", "make", "input"), ("make", "E", "output")])
    assert recycled_root_occurrences(["E"], {"E"}, inc, fwd) == []


def test_derived_in_edges_do_not_disqualify():
    inc, fwd = _graph([("C", "E", "dissociation"), ("K", "E", "depletion"), ("E", "R", "input")])
    assert recycled_root_occurrences(["E"], {"E"}, inc, fwd) == ["E"]


def test_a_dissociation_sink_is_never_a_root():
    # Review of PR #73: a sink fed only by dissociation, with no out-edges,
    # was pinned -- the case looked scored while the perturbation reached nothing.
    inc, fwd = _graph([("C", "SINK", "dissociation"), ("X", "C", "input")])
    assert recycled_root_occurrences(["SINK"], {"SINK"}, inc, fwd) == []


def test_complex_only_gene_pins_the_cycle_pool():
    # ALKBH2 shape: only complexes contain the gene, all in one regeneration cycle.
    inc, fwd = _graph([("EF", "bind", "input"), ("D", "bind", "input"), ("bind", "EFD", "output"),
                       ("EFD", "cat", "input"), ("cat", "EF", "output")])
    got = recycled_root_occurrences(["EF", "EFD"], set(), inc, fwd)
    assert sorted(got) == ["EF", "EFD"]


from benchmark_vs_mpbiopath import apply_catalog_ids  # noqa: E402


def test_catalog_list_overrides_the_mpbiopath_id(tmp_path):
    cat = tmp_path / "catalog_pathways.tsv"
    cat.write_text("# comment\nid\tpathway_name\nR-HSA-451927\tInterleukin-2_family_signaling\n"
                   "R-HSA-69620\tCell_Cycle_Checkpoints\n")
    got = apply_catalog_ids([("447115", "Interleukin-2_family_signaling"),
                             ("69620", "Cell_Cycle_Checkpoints"), ("1", "Not_listed")], cat)
    assert got == [("451927", "Interleukin-2_family_signaling"),
                   ("69620", "Cell_Cycle_Checkpoints"), ("1", "Not_listed")]


def test_pin_network_root_reads_pins_from_the_other_catalog(tmp_path, monkeypatch):
    import benchmark_vs_mpbiopath as bm
    real, served = tmp_path / "real" / "R-HSA-1", tmp_path / "served" / "R-HSA-1"
    for d, edges in ((real, "a,x\n"), (served, "x,a\n")):
        d.mkdir(parents=True)
        (d / "nodes.csv").write_text("uuid,node_kind\na,simple_entity\nx,simple_entity\n")
        (d / "logic_network.csv").write_text(
            "source_id,target_id,pos_neg,and_or,edge_type,stoichiometry,edge_reaction_id\n"
            + edges.replace("\n", ",pos,and,input,1,\n"))
    monkeypatch.setattr(bm, "PIN_SCOPE", "root")
    g2u = {"G": ["a", "x"]}
    monkeypatch.setattr(bm, "PIN_NETWORK_ROOT", "")
    assert bm.pin_occurrences(served, g2u, {})[0] == {"G": ["x"]}
    monkeypatch.setattr(bm, "PIN_NETWORK_ROOT", str(tmp_path / "real"))
    assert bm.pin_occurrences(served, g2u, {})[0] == {"G": ["a"]}
    monkeypatch.setattr(bm, "PIN_NETWORK_ROOT", str(tmp_path / "nowhere"))
    import pytest
    with pytest.raises(SystemExit):
        bm.pin_occurrences(served, g2u, {})


def test_pin_network_root_refuses_another_builds_nodes(tmp_path, monkeypatch):
    import pytest
    import benchmark_vs_mpbiopath as bm
    for name, ids in (("real", "a\nb\n"), ("served", "a\nc\n")):
        d = tmp_path / name / "R-HSA-1"
        d.mkdir(parents=True)
        (d / "logic_network.csv").write_text("source_id,target_id,pos_neg,and_or,edge_type,stoichiometry,edge_reaction_id\n")
        (d / "nodes.csv").write_text("uuid\n" + ids)
    monkeypatch.setattr(bm, "PIN_NETWORK_ROOT", str(tmp_path / "real"))
    with pytest.raises(SystemExit, match="does not list the served"):
        bm.pin_network_dir(tmp_path / "served" / "R-HSA-1")


def _depletion_net(tmp_path):
    # M (the gene's plain form) is fed only by a depletion edge from L; P (a
    # side form of the gene) has no in-edge at all. Both feed reactions.
    d = tmp_path / "R-HSA-1"
    d.mkdir()
    (d / "nodes.csv").write_text("uuid,node_kind\nm,simple_entity\np,simple_entity\nl,simple_entity\n"
                                 "r,reaction\nq,reaction\n")
    h = "source_id,target_id,pos_neg,and_or,edge_type,stoichiometry,edge_reaction_id\n"
    (d / "logic_network.csv").write_text(h + "l,m,neg,and,depletion,1,\nm,r,pos,and,input,1,R\n"
                                         "p,q,pos,and,input,1,Q\n")
    return d


@pytest.mark.parametrize("scope", ["root", "root_cycle"])
def test_pin_root_ignore_depletion(tmp_path, monkeypatch, scope):
    # specs/058: unset, M is not a root (its depletion in-edge counts) and only
    # P is pinned; under DS_PIN_ROOT_IGNORE=depletion both are.
    import benchmark_vs_mpbiopath as bm
    d = _depletion_net(tmp_path)
    monkeypatch.setattr(bm, "PIN_SCOPE", scope)
    monkeypatch.setattr(bm, "PIN_NETWORK_ROOT", "")
    monkeypatch.setattr(bm, "PIN_ROOT_IGNORE", "")
    assert bm.pin_occurrences(d, {"G": ["m", "p"]}, {})[0] == {"G": ["p"]}
    monkeypatch.setattr(bm, "PIN_ROOT_IGNORE", "depletion")
    assert bm.pin_occurrences(d, {"G": ["m", "p"]}, {})[0] == {"G": ["m", "p"]}


def test_pin_root_ignore_keeps_other_in_edges(tmp_path, monkeypatch):
    # A form with a real producer stays a non-root under the setting.
    import benchmark_vs_mpbiopath as bm
    d = _depletion_net(tmp_path)
    (d / "logic_network.csv").write_text((d / "logic_network.csv").read_text()
                                         + "q,p,pos,or,output,1,Q\n")
    monkeypatch.setattr(bm, "PIN_SCOPE", "root")
    monkeypatch.setattr(bm, "PIN_NETWORK_ROOT", "")
    monkeypatch.setattr(bm, "PIN_ROOT_IGNORE", "depletion")
    assert bm.pin_occurrences(d, {"G": ["m", "p"]}, {})[0] == {"G": ["m"]}


def test_pin_root_ignore_needs_an_out_edge(tmp_path, monkeypatch):
    # Week review 2026-10-15: a depletion-fed dead end (D: fed only by depletion,
    # feeds nothing) is not a root; pinning it would move nothing and block the
    # specs/024 fallback.
    import benchmark_vs_mpbiopath as bm
    d = _depletion_net(tmp_path)
    (d / "logic_network.csv").write_text((d / "logic_network.csv").read_text()
                                         + "l,dd,neg,and,depletion,1,\n")
    monkeypatch.setattr(bm, "PIN_SCOPE", "root")
    monkeypatch.setattr(bm, "PIN_NETWORK_ROOT", "")
    monkeypatch.setattr(bm, "PIN_ROOT_IGNORE", "depletion")
    assert bm.pin_occurrences(d, {"G": ["m", "dd"]}, {})[0] == {"G": ["m"]}


def test_pin_root_ignore_rejects_other_values():
    import subprocess
    r = subprocess.run([sys.executable, "-c", "import benchmark_vs_mpbiopath"],
                       cwd=str(Path(__file__).parent), capture_output=True, text=True,
                       env={**__import__("os").environ, "DS_PIN_ROOT_IGNORE": "dissociation"})
    assert r.returncode != 0 and "DS_PIN_ROOT_IGNORE" in r.stderr


def test_pin_root_ignore_default_is_depletion_and_none_restores_specs024():
    # specs/061 adoption (2026-10-15): depletion-only forms are roots by default
    import os
    import subprocess
    env = {k: v for k, v in os.environ.items() if k != "DS_PIN_ROOT_IGNORE"}
    probe = "import benchmark_vs_mpbiopath as bm; print(repr(bm.PIN_ROOT_IGNORE))"
    for value, want in ((None, "'depletion'"), ("none", "''"), ("depletion", "'depletion'")):
        e = env if value is None else {**env, "DS_PIN_ROOT_IGNORE": value}
        r = subprocess.run([sys.executable, "-c", probe], cwd=str(Path(__file__).parent),
                           capture_output=True, text=True, env=e)
        assert r.returncode == 0, r.stderr
        assert r.stdout.strip().splitlines()[-1] == want


def test_a_recycle_in_edge_does_not_stop_a_root(tmp_path, monkeypatch):
    # specs/081: free p130 (f) -> starting copy (s) by a recycle edge; s builds the
    # complex c. s is still the root a p130 knockout pins, under root and root_cycle.
    import benchmark_vs_mpbiopath as bm
    d = tmp_path / "R-HSA-1"
    d.mkdir()
    (d / "nodes.csv").write_text("uuid,node_kind,diagram_entity_id\nf,simple_entity,R-P\n"
                                 "s,simple_entity,R-P\nc,simple_complex,R-C\nr,reaction,R-R\n")
    rows = [("s", "c", "assembly"), ("c", "r", "input"), ("r", "f", "output"), ("f", "s", "recycle")]
    (d / "logic_network.csv").write_text(
        "source_id,target_id,pos_neg,and_or,edge_type,stoichiometry,edge_reaction_id\n"
        + "".join(f"{a},{b},pos,and,{t},1,\n" for a, b, t in rows))
    monkeypatch.setattr(bm, "PIN_NETWORK_ROOT", "")
    for scope in ("root", "root_cycle"):
        monkeypatch.setattr(bm, "PIN_SCOPE", scope)
        assert bm.pin_occurrences(d, {"P": ["f", "s"]}, {})[0] == {"P": ["s"]}
