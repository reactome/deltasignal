"""specs/051 protocol: perturbation-clustered bootstrap intervals and grouped folds.

Cases that share a perturbation are correlated (one drug solve reads many sites),
so intervals resample whole groups, never cases. Folds are assigned per group
from a hash of (seed, group), so they are stable across runs, machines and
case order, and committing the seed fixes them.
"""
import hashlib
import random


def cluster_bootstrap_ci(rows, stat, cluster, n_boot=2000, seed=51, alpha=0.05):
    """(point, low, high) for stat(rows), resampling whole clusters with replacement.

    `stat` maps a list of rows to a number; `cluster` maps a row to its group key."""
    groups = {}
    for r in rows:
        groups.setdefault(cluster(r), []).append(r)
    keys = sorted(groups, key=str)
    point = stat(rows)
    if len(keys) < 2:
        # one cluster carries no between-cluster variation: no interval, not a
        # zero-width one that a "low > 0" test would pass (code review 2026-10-13)
        return point, None, None
    rng = random.Random(seed)
    draws = []
    for _ in range(n_boot):
        sample = [r for _k in range(len(keys)) for r in groups[keys[rng.randrange(len(keys))]]]
        draws.append(stat(sample))
    draws.sort()
    # symmetric order statistics: as many draws below `lo` as above `hi`
    lo = draws[round(alpha / 2 * (n_boot - 1))]
    hi = draws[round((1 - alpha / 2) * (n_boot - 1))]
    return point, lo, hi


def assign_fold(group, k, seed):
    """Fold 0..k-1 for a group key, from sha256(seed:group): order- and run-independent."""
    h = hashlib.sha256(f"{seed}:{group}".encode()).hexdigest()
    return int(h, 16) % k


def folds(groups, k, seed):
    """{group: fold} for every distinct group."""
    return {g: assign_fold(g, k, seed) for g in sorted(set(groups), key=str)}
