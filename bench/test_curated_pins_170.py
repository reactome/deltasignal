"""specs/170: under DS_PIN_CURATED_ROOT a literature edge cannot stop a curated gene being pinned; literature nodes are
rooted by the network as served. Runs the real pin_occurrences on fixture pathway directories."""
import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
import benchmark_vs_mpbiopath as bm  # noqa: E402

FIELDS = ["source_id", "target_id", "pos_neg", "and_or", "edge_type", "stoichiometry", "edge_reaction_id"]
P = "R-HSA-0000001"


def _write(d, edges, nodes):
    d.mkdir(parents=True)
    with open(d / "logic_network.csv", "w", newline="") as f:
        w = csv.DictWriter(f, FIELDS)
        w.writeheader()
        for s, t, et in edges:
            w.writerow({"source_id": s, "target_id": t, "pos_neg": "pos", "and_or": "and", "edge_type": et,
                        "stoichiometry": "1", "edge_reaction_id": ""})
    with open(d / "nodes.csv", "w", newline="") as f:
        w = csv.DictWriter(f, ["uuid", "diagram_entity_id"])
        w.writeheader()
        for u in nodes:
            w.writerow({"uuid": u, "diagram_entity_id": u})


# Curated: ERBB2 -> R1 -> PERBB2 -> R2 -> AKT_OUT. Literature: PHF8 -> LR -> ERBB2 (into a curated root),
# and LIT2 fed by LR (a literature node with an in-edge, so not a root), plus PERBB2 -> LIT3 (out of a curated node).
CURATED = [("ERBB2", "R1", "input"), ("R1", "PERBB2", "output"), ("PERBB2", "R2", "input"),
           ("R2", "AKT_OUT", "output")]
LITERATURE = [("PHF8", "LR", "input"), ("LR", "ERBB2", "output"), ("LR", "LIT2", "output"),
              ("PERBB2", "LIT3", "input")]
CURATED_NODES = ["ERBB2", "R1", "PERBB2", "R2", "AKT_OUT"]


@pytest.fixture
def arms(tmp_path, monkeypatch):
    served = tmp_path / "with" / P
    curated = tmp_path / "without" / P
    _write(served, CURATED + LITERATURE, CURATED_NODES + ["PHF8", "LR", "LIT2", "LIT3"])
    _write(curated, CURATED, CURATED_NODES)
    monkeypatch.setattr(bm, "PIN_SCOPE", "root_cycle")
    monkeypatch.setattr(bm, "PIN_NETWORK_ROOT", "")
    return served, tmp_path / "without"


def _pins(served, genes):
    pins, _ = bm.pin_occurrences(served, {g: [g] for g in genes}, {g: [g] for g in genes})
    return pins


def test_without_the_flag_a_literature_edge_unroots_the_curated_gene(arms, monkeypatch):
    served, _ = arms
    monkeypatch.setattr(bm, "PIN_CURATED_ROOT", "")
    # ERBB2 is fed by LR; no own-cycle fallback applies, so it is not perturbed (the specs/168 mechanism)
    assert _pins(served, ["ERBB2"])["ERBB2"] == []


def test_with_the_flag_the_curated_gene_is_pinned_as_the_curated_network_pins_it(arms, monkeypatch):
    served, curated_root = arms
    monkeypatch.setattr(bm, "PIN_CURATED_ROOT", str(curated_root))
    assert _pins(served, ["ERBB2"])["ERBB2"] == ["ERBB2"]


def test_literature_nodes_are_rooted_by_the_combined_network(arms, monkeypatch):
    served, curated_root = arms
    monkeypatch.setattr(bm, "PIN_CURATED_ROOT", str(curated_root))
    pins = _pins(served, ["PHF8", "LIT2"])
    assert pins["PHF8"] == ["PHF8"]      # fed by nothing: a root
    assert pins["LIT2"] == []            # fed by a literature reaction: not a root


def test_a_curated_mid_pathway_node_stays_unpinned(arms, monkeypatch):
    served, curated_root = arms
    monkeypatch.setattr(bm, "PIN_CURATED_ROOT", str(curated_root))
    assert _pins(served, ["PERBB2"])["PERBB2"] == []


def test_the_filter_drops_only_edges_into_curated_nodes_the_curated_network_lacks(arms, monkeypatch):
    served, curated_root = arms
    monkeypatch.setattr(bm, "PIN_CURATED_ROOT", str(curated_root))
    with open(served / "logic_network.csv") as f:
        kept = {(r["source_id"], r["target_id"]) for r in bm.curated_pin_rows(served, csv.DictReader(f))}
    assert ("LR", "ERBB2") not in kept                     # literature edge into a curated node
    assert {("PHF8", "LR"), ("LR", "LIT2"), ("PERBB2", "LIT3")} <= kept
    assert {(s, t) for s, t, _ in CURATED} <= kept


def test_a_curated_network_with_foreign_nodes_is_refused(arms, monkeypatch, tmp_path):
    served, _ = arms
    other = tmp_path / "other"
    _write(other / P, CURATED, CURATED_NODES + ["NOT_SERVED"])
    monkeypatch.setattr(bm, "PIN_CURATED_ROOT", str(other))
    with pytest.raises(SystemExit, match="lacks"):
        _pins(served, ["ERBB2"])


def test_a_missing_curated_pathway_is_refused(arms, monkeypatch, tmp_path):
    served, _ = arms
    monkeypatch.setattr(bm, "PIN_CURATED_ROOT", str(tmp_path / "nowhere"))
    with pytest.raises(SystemExit, match="no R-HSA-0000001"):
        _pins(served, ["ERBB2"])
