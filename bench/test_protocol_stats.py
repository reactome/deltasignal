"""specs/051: clustered bootstrap and grouped folds."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import protocol_stats as ps  # noqa: E402


def acc(rows):
    return sum(r["ok"] for r in rows) / len(rows) if rows else 0.0


def test_bootstrap_resamples_clusters_not_cases():
    # one drug with 100 identical right cases, one drug with 100 wrong: by cluster,
    # accuracy is either 0, 0.5 or 1 in every draw; by case it would sit near 0.5
    rows = [{"d": "a", "ok": 1}] * 100 + [{"d": "b", "ok": 0}] * 100
    point, lo, hi = ps.cluster_bootstrap_ci(rows, acc, lambda r: r["d"], n_boot=500)
    assert point == 0.5 and lo == 0.0 and hi == 1.0


def test_bootstrap_deterministic_and_single_cluster():
    rows = [{"d": i % 7, "ok": i % 3 == 0} for i in range(70)]
    a = ps.cluster_bootstrap_ci(rows, acc, lambda r: r["d"], n_boot=300, seed=3)
    b = ps.cluster_bootstrap_ci(list(reversed(rows)), acc, lambda r: r["d"], n_boot=300, seed=3)
    assert a == b                                   # case order does not matter
    assert a[1] <= a[0] <= a[2]
    one = [{"d": "x", "ok": 1}, {"d": "x", "ok": 0}]
    assert ps.cluster_bootstrap_ci(one, acc, lambda r: r["d"]) == (0.5, 0.5, 0.5)


def test_folds_stable_and_grouped():
    f1 = ps.folds(["EGFR", "SRC", "EGFR", "MTOR"], 5, seed=51)
    f2 = ps.folds(["MTOR", "SRC", "EGFR"], 5, seed=51)
    assert f1 == f2 and set(f1) == {"EGFR", "SRC", "MTOR"}
    assert all(0 <= v < 5 for v in f1.values())
    assert ps.assign_fold("EGFR", 5, 51) == f1["EGFR"]
    spread = ps.folds([f"g{i}" for i in range(1000)], 5, seed=51)
    counts = [list(spread.values()).count(i) for i in range(5)]
    assert min(counts) > 150                        # roughly balanced
