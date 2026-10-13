"""specs/051: the rewired-network null keeps degrees, classes and node identity."""
import csv
import random
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rewire_network as rw  # noqa: E402


def _row(s, t, et="input", pn="pos", ao="and", rid="R1"):
    return {"source_id": s, "target_id": t, "pos_neg": pn, "and_or": ao,
            "edge_type": et, "stoichiometry": "1", "edge_reaction_id": rid}


def _net(n=40, seed=0):
    rng = random.Random(seed)
    ents = [f"e{i}" for i in range(n)]
    rxns = [f"r{i}" for i in range(n)]
    rows, seen = [], set()
    while len(rows) < 3 * n:
        s, t = rng.choice(ents), rng.choice(rxns)
        if (s, t) not in seen:
            seen.add((s, t))
            rows.append(_row(s, t))
    for i in range(n):
        rows.append(_row(rxns[i], ents[(i + 1) % n], et="output", ao="or"))
        if i % 3 == 0:
            rows.append(_row(ents[i], rxns[(i + 5) % n], et="regulator", pn="neg", ao="or"))
    rows.append(_row("e0", "v0", et="variant_pool", ao="or"))
    return rows


def test_swap_class_preserves_degrees_and_rejects_duplicates():
    pairs = [("a", "x"), ("b", "y"), ("c", "z"), ("a", "y")]
    before = list(pairs)
    acc, tried = rw.swap_class(pairs, random.Random(1), swaps_per_edge=5, tries_per_edge=50)
    assert Counter(s for s, _ in pairs) == Counter(s for s, _ in before)
    assert Counter(t for _, t in pairs) == Counter(t for _, t in before)
    assert len(set(pairs)) == len(pairs)          # no repeated pair
    assert all(s != t for s, t in pairs)
    assert 0 < acc <= tried
    assert rw.swap_class([("a", "x")], random.Random(1)) == (0, 0)


def test_rewire_rows_keeps_classes_and_degrees_and_moves_edges():
    rows = _net()
    new, stats = rw.rewire_rows(rows, seed=1, pathway="P")
    assert rw.degrees(new) == rw.degrees(rows)
    assert [rw.edge_class(r) for r in new] == [rw.edge_class(r) for r in rows]
    assert [r["source_id"] for r in new] == [r["source_id"] for r in rows]
    # an entity-structure edge is never touched
    assert new[-1] == rows[-1]
    moved = sum(r["target_id"] != o["target_id"] for r, o in zip(new, rows))
    assert moved == sum(s["moved"] for s in stats.values()) > len(rows) // 3
    # a moved edge loses its stale reaction id; an unmoved one keeps it
    for r, o in zip(new, rows):
        assert r["edge_reaction_id"] == ("" if r["target_id"] != o["target_id"] else o["edge_reaction_id"])
    assert set(stats) == {"input|pos|and", "output|pos|or", "regulator|neg|or"}


def test_rewire_is_reproducible_by_seed_and_varies_across_seeds():
    rows = _net()
    a, _ = rw.rewire_rows(rows, 1, "P")
    b, _ = rw.rewire_rows(list(reversed(rows)), 1, "P")
    key = lambda rs: sorted((r["source_id"], r["target_id"], r["edge_type"]) for r in rs)  # noqa: E731
    assert key(a) == key(b)                       # input order does not matter
    c, _ = rw.rewire_rows(rows, 2, "P")
    assert key(a) != key(c)


def test_rewire_pathway_drops_wiring_tables(tmp_path):
    src, dst = tmp_path / "R-HSA-1", tmp_path / "out"
    src.mkdir()
    (src / "cache").mkdir()
    fields = list(_row("a", "b"))
    rw.write_network(src / "logic_network.csv", fields, _net(10))
    for name in rw.DROPPED_TABLES + ("containment.csv", "nodes.csv"):
        (src / name).write_text("h\n")
    dst.mkdir()
    rw.rewire_pathway(src, dst / "R-HSA-1", seed=3)
    names = {p.name for p in (dst / "R-HSA-1").iterdir()}
    assert names == {"logic_network.csv", "containment.csv", "nodes.csv"}
    with open(dst / "R-HSA-1" / "logic_network.csv") as f:
        assert csv.DictReader(f).fieldnames == fields
