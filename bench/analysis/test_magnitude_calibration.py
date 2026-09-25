"""Tests for magnitude_calibration.py.

The statistic decides whether a claim goes in the paper, so these pin the two
ways it could lie: an AUC that is not an AUC, and stratification that does not
actually remove composition.
"""
import os
import sys

import math

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


# --- Pins added after review: eight mutants survived the first version. ---

from magnitude_calibration import (  # noqa: E402
    load, report, shuffle_within, split_rows,
)
from holdout_report import TUNING_PATHWAYS  # noqa: E402


def test_shuffle_is_within_strata_not_across():
    """Each stratum must keep exactly its own values. Shuffling across strata
    would destroy the composition the null is meant to preserve."""
    strata = [(np.array([1.0, 2.0, 3.0]), np.array([True, False, True])),
              (np.array([10.0, 20.0]), np.array([False, True]))]
    out = shuffle_within(strata, np.random.default_rng(0))
    for (s0, c0), (s1, c1) in zip(strata, out):
        assert sorted(s0) == sorted(s1)                       # same values
        assert (c0 == c1).all()                               # outcomes untouched


def test_p_is_never_zero():
    signal = build_strata(
        [{"pathway": "A", "cls": "0", "correct": i >= 20, "s": float(i)} for i in range(40)],
        lambda r: r["pathway"])
    obs, _ = stratified_auc(signal)
    assert permutation_p(signal, obs, 10, 0) >= 1 / 11


def test_ties_with_the_observed_count_as_at_least_as_good():
    """All strengths tied: every permutation equals the observed 0.5, so with
    `>=` the p is exactly 1. A `>` would report 1/(n+1) -- a false signal."""
    tied = build_strata(
        [{"pathway": "A", "cls": "0", "correct": i % 2 == 0, "s": 1.0} for i in range(20)],
        lambda r: r["pathway"])
    obs, _ = stratified_auc(tied)
    assert obs == 0.5
    assert permutation_p(tied, obs, 50, 0) == 1.0


def test_strata_are_combined_by_pairs_not_by_equal_weight():
    """The pre-registered statistic weights strata by correct x incorrect pairs.
    This pins that choice: a big flat stratum must dominate a small perfect one."""
    big_flat = (np.array([1.0] * 20 + [1.0] * 20), np.array([True] * 20 + [False] * 20))
    small_perfect = (np.array([2.0, 1.0]), np.array([True, False]))
    auc, pairs = stratified_auc([big_flat, small_perfect])
    assert pairs == 400 + 1
    assert abs(auc - (0.5 * 400 + 1.0 * 1) / 401) < 1e-12     # pair-weighted
    assert abs(auc - 0.75) > 0.2                              # not the equal-weight mean


def _write(path, rows):
    cols = ["pathway", "gene", "direction", "key_output", "predicted", "expected", "pred_ui"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")


def test_load_excludes_normal_and_scores_correctness(tmp_path):
    p = tmp_path / "cases.tsv"
    _write(p, [
        {"pathway": "P", "gene": "G", "direction": "0", "key_output": "1",
         "predicted": "0", "expected": "0", "pred_ui": "0.1"},       # DOWN, right
        {"pathway": "P", "gene": "G", "direction": "0", "key_output": "2",
         "predicted": "2", "expected": "0", "pred_ui": "5"},         # UP, wrong
        {"pathway": "P", "gene": "G", "direction": "0", "key_output": "3",
         "predicted": "1", "expected": "1", "pred_ui": "1"},         # NORMAL: excluded
    ])
    rows = load(str(p))
    assert len(rows) == 2
    assert [r["correct"] for r in rows] == [True, False]
    assert abs(rows[0]["s"] - 1.0) < 1e-12 and abs(rows[1]["s"] - math.log10(5)) < 1e-12


def test_split_puts_tuning_pathways_on_the_tuning_side():
    tuning_name = sorted(TUNING_PATHWAYS)[0]
    rows = [{"pathway": tuning_name}, {"pathway": "Some_held_out_pathway"}]
    held, tune = split_rows(rows)
    assert [r["pathway"] for r in tune] == [tuning_name]
    assert [r["pathway"] for r in held] == ["Some_held_out_pathway"]


def test_report_decides_on_pathway_and_class_strata():
    """The decision is made on (pathway, class) strata. One pathway with two
    predicted classes must yield two strata; a pathway-only stratification
    would pool them into one."""
    rows = []
    for i in range(20):
        rows.append({"pathway": "P", "cls": "0", "correct": i >= 10, "s": float(i)})        # signal
        rows.append({"pathway": "P", "cls": "2", "correct": i >= 10, "s": float(i)})        # signal
    r = report("x", rows, 20, 0)
    assert r["within (pathway, class)"]["strata"] == 2
    assert r["within (pathway, class)"]["auc"] > 0.95
