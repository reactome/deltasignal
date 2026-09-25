"""Tests for magnitude_calibration.py.

The statistic decides whether a claim goes in the paper, so these pin the two
ways it could lie: an AUC that is not an AUC, and stratification that does not
actually remove composition.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from magnitude_calibration import (  # noqa: E402
    build_strata,
    mann_whitney_parts,
    permutation_p,
    strength,
    stratified_auc,
)


def test_strength_is_symmetric_in_log_space_and_floors_zero():
    assert strength(1.0) == 0.0
    assert abs(strength(10.0) - strength(0.1)) < 1e-12      # 10x up == 10x down
    assert strength(100.0) == 2.0
    assert strength(0.0) == 6.0                             # floored, not undefined


def test_auc_is_one_zero_and_half_at_the_extremes():
    c = np.array([True, True, False, False])
    u, pairs = mann_whitney_parts(np.array([3.0, 4.0, 1.0, 2.0]), c)
    assert (u / pairs) == 1.0                               # correct all stronger
    u, pairs = mann_whitney_parts(np.array([1.0, 2.0, 3.0, 4.0]), c)
    assert (u / pairs) == 0.0                               # correct all weaker
    u, pairs = mann_whitney_parts(np.array([5.0, 5.0, 5.0, 5.0]), c)
    assert (u / pairs) == 0.5                               # all tied


def test_auc_matches_a_brute_force_pair_count():
    rng = np.random.default_rng(1)
    s = rng.integers(0, 5, 40).astype(float)                # plenty of ties
    c = rng.random(40) < 0.5
    u, pairs = mann_whitney_parts(s, c)
    brute = sum((a > b) + 0.5 * (a == b) for a in s[c] for b in s[~c])
    assert abs(u - brute) < 1e-9
    assert pairs == c.sum() * (~c).sum()


def test_strata_with_one_outcome_are_uninformative_and_dropped():
    rows = [{"pathway": "P", "cls": "0", "correct": True, "s": 1.0},
            {"pathway": "P", "cls": "0", "correct": True, "s": 2.0}]
    assert build_strata(rows, lambda r: (r["pathway"], r["cls"])) == []


def test_stratification_removes_pure_composition():
    """The whole point. Pathway A is always right and always strong; pathway B
    always wrong and always weak -- but WITHIN each pathway strength says
    nothing. Pooled AUC is a perfect 1.0; stratified must be 0.5."""
    rows = []
    rng = np.random.default_rng(7)
    for i in range(60):
        rows.append({"pathway": "A", "cls": "0", "correct": i % 4 != 0, "s": 5 + rng.random()})
        rows.append({"pathway": "B", "cls": "0", "correct": i % 4 == 0, "s": 1 + rng.random()})
    pooled, _ = stratified_auc(build_strata(rows, lambda r: 0))
    within, _ = stratified_auc(build_strata(rows, lambda r: r["pathway"]))
    assert pooled > 0.7                                     # composition looks like signal
    assert abs(within - 0.5) < 0.12                         # and stratifying removes it


def test_real_within_stratum_signal_survives_stratification():
    rows = []
    for p in ("A", "B", "C"):
        for i in range(40):
            rows.append({"pathway": p, "cls": "0", "correct": i >= 20, "s": float(i)})
    within, _ = stratified_auc(build_strata(rows, lambda r: r["pathway"]))
    assert within > 0.95


def test_permutation_p_is_small_for_signal_and_large_for_noise():
    signal = build_strata(
        [{"pathway": "A", "cls": "0", "correct": i >= 20, "s": float(i)} for i in range(40)],
        lambda r: r["pathway"])
    obs, _ = stratified_auc(signal)
    assert permutation_p(signal, obs, 300, 0) < 0.01
    rng = np.random.default_rng(3)
    noise = build_strata(
        [{"pathway": "A", "cls": "0", "correct": bool(rng.random() < 0.5), "s": float(rng.random())}
         for _ in range(40)],
        lambda r: r["pathway"])
    obs, _ = stratified_auc(noise)
    assert permutation_p(noise, obs, 300, 0) > 0.05
