#!/usr/bin/env python3
"""specs/127: the registered scoring of Lun et al. 2019, as pure functions.

Each function is one numbered step of specs/127-lun2019-scoring/spec.md, tested on synthetic cells in
bench/test_lun_score.py. Nothing here opens the dataset. The loader that turns `complete_dat.csv` into the inputs below
is written after the format check (which reads the header and sample identifiers only) and committed before any marker
value is read.

Inputs, per sample: GFP per cell and readout values per cell (lists of floats, as provided; step 1 transforms them).
"""
import math
import random
from collections import defaultdict

import numpy as np

UP, DOWN, NOT, CONFLICT = "up", "down", "not", "conflict"
CLASSES = (UP, DOWN, NOT)
COFACTOR = 5.0             # step 1: arcsinh(x / 5)
TRANSFECTED_PCT = 99.0     # step 2
N_BINS = 5                 # step 3: GFP quintiles of transfected FLAG-GFP cells
MIN_BIN_CELLS = 20         # step 4
MIN_BINS = 2
MIN_SAMPLE_CELLS = 100
NULL_PCT = 95.0            # step 5
MIN_CONTROLS = 10
BOOTSTRAP = 10_000


def asinh(values):
    """Step 1."""
    return np.arcsinh(np.asarray(values, dtype=float) / COFACTOR)


def transfected_cutoff(untransfected_gfp):
    """Step 2: a cell is transfected if its (transformed) GFP exceeds this."""
    return float(np.percentile(asinh(untransfected_gfp), TRANSFECTED_PCT))


def bin_edges(control_gfp, cutoff):
    """Step 3: the interior quintile edges of transfected FLAG-GFP cells' (transformed) GFP."""
    g = asinh(control_gfp)
    g = g[g > cutoff]
    if g.size < N_BINS:
        raise ValueError(f"only {g.size} transfected FLAG-GFP cells: cannot form {N_BINS} dose bins")
    return [float(x) for x in np.percentile(g, [100.0 * k / N_BINS for k in range(1, N_BINS)])]


def _binned(gfp, values, cutoff, edges):
    """{bin: transformed readout values} over transfected cells."""
    g, v = asinh(gfp), asinh(values)
    keep = g > cutoff
    bins = np.searchsorted(edges, g[keep], side="right")
    out = defaultdict(list)
    for b, x in zip(bins, v[keep]):
        out[int(b)].append(float(x))
    return out


def n_transfected(gfp, cutoff):
    return int((asinh(gfp) > cutoff).sum())


def effect(sample_gfp, sample_values, control_gfp, control_values, cutoff, edges):
    """Step 4: the median over dose bins of (sample median − control median); None if fewer than MIN_BINS bins have
    MIN_BIN_CELLS sample cells (the control side needs at least one cell in the bin)."""
    s = _binned(sample_gfp, sample_values, cutoff, edges)
    c = _binned(control_gfp, control_values, cutoff, edges)
    d = [float(np.median(s[b])) - float(np.median(c[b]))
         for b in sorted(s) if len(s[b]) >= MIN_BIN_CELLS and c.get(b)]
    return float(np.median(d)) if len(d) >= MIN_BINS else None


def protein_effect(effects):
    """Step 4: a protein's E is the median over its samples' E (None ignored); None if no sample has one."""
    e = [x for x in effects if x is not None]
    return float(np.median(e)) if e else None


def null_threshold(controls, cutoff, edges):
    """Step 5. `controls`: [(gfp, values)] for one readout, one entry per FLAG-GFP control sample. Each sample's E is
    taken against the pooled others (leave one out); δ is the 95th percentile of |E|, or the largest |E| when there
    are fewer than MIN_CONTROLS samples with an E. None if no control sample has an E."""
    es = []
    for i, (g, v) in enumerate(controls):
        og = [x for j, (gg, _) in enumerate(controls) if j != i for x in gg]
        ov = [x for j, (_, vv) in enumerate(controls) if j != i for x in vv]
        e = effect(g, v, og, ov, cutoff, edges)
        if e is not None:
            es.append(abs(e))
    if not es:
        return None
    return float(np.percentile(es, NULL_PCT)) if len(es) >= MIN_CONTROLS else max(es)


def measured_call(e, delta):
    """Step 6."""
    if e is None or delta is None:
        return None
    if e > delta:
        return UP
    if e < -delta:
        return DOWN
    return NOT


def positive_control(self_calls):
    """Step 7: (passes, n up, n measured): at least half of the measured self pairs read UP."""
    measured = [c for c in self_calls if c is not None]
    ups = sum(1 for c in measured if c == UP)
    return (bool(measured) and 2 * ups >= len(measured)), ups, len(measured)


def combine_networks(calls):
    """The model's prediction for a pair over the networks that score it: the call every change-predicting network
    agrees on; CONFLICT if they disagree; NOT if none predicts a change."""
    changed = {c for c in calls if c in (UP, DOWN)}
    if not changed:
        return NOT
    return changed.pop() if len(changed) == 1 else CONFLICT


def as_scored(pred):
    """A missing prediction or a conflict counts as NOT."""
    return pred if pred in (UP, DOWN) else NOT


def balanced_accuracy(pairs):
    """Three-class balanced accuracy: mean over measured classes present of the share predicted correctly.
    `pairs`: [(measured, predicted)]."""
    hit, tot = defaultdict(int), defaultdict(int)
    for m, p in pairs:
        tot[m] += 1
        hit[m] += int(as_scored(p) == m)
    present = [c for c in CLASSES if tot[c]]
    return sum(hit[c] / tot[c] for c in present) / len(present) if present else float("nan")


def bootstrap_difference(cases, model_a, model_b, n=BOOTSTRAP, seed=127):
    """Balanced accuracy of A minus B, with a 95% percentile interval from resampling proteins.
    `cases`: [{"gene", "measured", model_a: call, model_b: call}]. Returns (point, low, high)."""
    by = defaultdict(list)
    for c in cases:
        by[c["gene"]].append(c)
    genes = sorted(by)

    def diff(sample):
        return (balanced_accuracy([(c["measured"], c[model_a]) for c in sample])
                - balanced_accuracy([(c["measured"], c[model_b]) for c in sample]))

    point = diff(cases)
    rng = random.Random(seed)
    draws = []
    for _ in range(n):
        sample = [c for g in (rng.choice(genes) for _ in genes) for c in by[g]]
        d = diff(sample)
        if not math.isnan(d):
            draws.append(d)
    low, high = np.percentile(draws, [2.5, 97.5])
    return point, float(low), float(high)


def verdict(low, high):
    """Registered reading of an interval."""
    return "gain" if low > 0 else "loss" if high < 0 else "unresolved"
