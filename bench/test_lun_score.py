"""specs/127: the registered Lun 2019 scoring steps on synthetic cells."""
import math
import random
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lun_score as ls  # noqa: E402


def raw(z):
    """Values whose arcsinh(x / 5) is z, so shifts can be set on the transformed scale."""
    return list(5.0 * np.sinh(np.asarray(z, dtype=float)))


def cells(n, gfp_mu, readout_mu, seed):
    rng = np.random.default_rng(seed)
    return raw(rng.normal(gfp_mu, 0.6, n)), raw(rng.normal(readout_mu, 0.3, n))


@pytest.fixture
def setup():
    untrans_gfp, _ = cells(2000, 0.0, 2.0, 1)
    cutoff = ls.transfected_cutoff(untrans_gfp)
    ctrl_gfp, ctrl_val = cells(5000, 4.0, 2.0, 2)
    edges = ls.bin_edges(ctrl_gfp, cutoff)
    return cutoff, edges, ctrl_gfp, ctrl_val


def test_asinh_cofactor_five():
    assert ls.asinh([0.0, 5.0])[1] == pytest.approx(math.asinh(1.0))
    assert ls.asinh([-5.0])[0] == pytest.approx(-math.asinh(1.0))


def test_transfected_cutoff_is_99th_percentile_of_untransfected():
    z = np.linspace(-1, 1, 1001)
    assert ls.transfected_cutoff(raw(z)) == pytest.approx(np.percentile(z, 99))


def test_bin_edges_are_quintiles_of_transfected_controls(setup):
    cutoff, edges, ctrl_gfp, _ = setup
    assert len(edges) == ls.N_BINS - 1 and edges == sorted(edges) and edges[0] > cutoff
    g = ls.asinh(ctrl_gfp)
    g = g[g > cutoff]
    shares = np.histogram(g, bins=[-np.inf] + edges + [np.inf])[0] / g.size
    assert shares == pytest.approx([0.2] * 5, abs=0.01)


def test_bin_edges_refuse_too_few_transfected_cells():
    with pytest.raises(ValueError):
        ls.bin_edges(raw([5.0, 5.0]), 4.0)


def test_effect_recovers_a_shift_on_the_transformed_scale(setup):
    cutoff, edges, cg, cv = setup
    sg, sv = cells(3000, 4.0, 2.5, 3)
    assert ls.effect(sg, sv, cg, cv, cutoff, edges) == pytest.approx(0.5, abs=0.05)
    sg, sv = cells(3000, 4.0, 1.6, 4)
    assert ls.effect(sg, sv, cg, cv, cutoff, edges) == pytest.approx(-0.4, abs=0.05)


def test_effect_is_dose_matched(setup):
    """A sample whose cells sit at higher GFP, with a readout that tracks GFP in the controls too, has no effect."""
    cutoff, edges, _, _ = setup
    rng = np.random.default_rng(5)
    cz = rng.normal(4.0, 0.6, 5000)
    cg, cv = raw(cz), raw(1.0 + 0.5 * cz + rng.normal(0, 0.1, 5000))
    sz = rng.normal(4.6, 0.6, 3000)
    sg, sv = raw(sz), raw(1.0 + 0.5 * sz + rng.normal(0, 0.1, 3000))
    unmatched = float(np.median(ls.asinh(sv)) - np.median(ls.asinh(cv)))
    assert unmatched > 0.2
    assert abs(ls.effect(sg, sv, cg, cv, cutoff, edges)) < 0.1


def test_effect_needs_two_bins_of_twenty_cells(setup):
    cutoff, edges, cg, cv = setup
    lo = (edges[0] + cutoff) / 2          # all in bin 0
    sg, sv = raw([lo] * 50), raw([3.0] * 50)
    assert ls.effect(sg, sv, cg, cv, cutoff, edges) is None
    mid = (edges[1] + edges[2]) / 2
    sg2, sv2 = raw([lo] * 20 + [mid] * 20), raw([3.0] * 40)
    assert ls.effect(sg2, sv2, cg, cv, cutoff, edges) is not None
    sg3, sv3 = raw([lo] * 20 + [mid] * 19), raw([3.0] * 39)
    assert ls.effect(sg3, sv3, cg, cv, cutoff, edges) is None


def test_untransfected_cells_are_ignored(setup):
    cutoff, edges, cg, cv = setup
    sg, sv = cells(3000, 4.0, 2.0, 6)
    ug, uv = raw([0.0] * 3000), raw([6.0] * 3000)
    base = ls.effect(sg, sv, cg, cv, cutoff, edges)
    assert ls.effect(sg + ug, sv + uv, cg, cv, cutoff, edges) == pytest.approx(base)


def test_protein_effect_is_median_over_samples():
    assert ls.protein_effect([0.1, None, 0.5, 0.3]) == pytest.approx(0.3)
    assert ls.protein_effect([None]) is None


def test_null_threshold_from_controls_separates_a_real_effect(setup):
    cutoff, edges, _, _ = setup
    controls = [cells(1500, 4.0, 2.0, 100 + i) for i in range(12)]
    delta = ls.null_threshold(controls, cutoff, edges)
    assert 0 < delta < 0.1
    sg, sv = cells(1500, 4.0, 2.4, 7)
    pooled_g = [x for g, _ in controls for x in g]
    pooled_v = [x for _, v in controls for x in v]
    e = ls.effect(sg, sv, pooled_g, pooled_v, cutoff, edges)
    assert ls.measured_call(e, delta) == ls.UP
    sg, sv = cells(1500, 4.0, 2.0, 8)
    e0 = ls.effect(sg, sv, pooled_g, pooled_v, cutoff, edges)
    assert ls.measured_call(e0, delta) == ls.NOT


def test_null_threshold_with_few_controls_is_the_largest(setup, monkeypatch):
    cutoff, edges, _, _ = setup
    controls = [cells(1500, 4.0, 2.0, 200 + i) for i in range(4)]
    monkeypatch.setattr(ls, "NULL_PCT", 50.0)
    small = ls.null_threshold(controls, cutoff, edges)
    es = []
    for i, (g, v) in enumerate(controls):
        og = [x for j, (gg, _) in enumerate(controls) if j != i for x in gg]
        ov = [x for j, (_, vv) in enumerate(controls) if j != i for x in vv]
        es.append(abs(ls.effect(g, v, og, ov, cutoff, edges)))
    assert small == pytest.approx(max(es))


def test_measured_call_is_strict():
    assert ls.measured_call(0.2, 0.2) == ls.NOT
    assert ls.measured_call(0.21, 0.2) == ls.UP
    assert ls.measured_call(-0.21, 0.2) == ls.DOWN
    assert ls.measured_call(None, 0.2) is None
    assert ls.measured_call(0.5, None) is None


def test_positive_control_half_up():
    assert ls.positive_control([ls.UP, ls.UP, ls.NOT, ls.DOWN, None]) == (True, 2, 4)
    assert ls.positive_control([ls.UP, ls.NOT, ls.NOT]) == (False, 1, 3)
    assert ls.positive_control([None]) == (False, 0, 0)


def test_combine_networks():
    assert ls.combine_networks([ls.NOT, ls.UP, ls.UP]) == ls.UP
    assert ls.combine_networks([ls.DOWN, ls.UP]) == ls.CONFLICT
    assert ls.combine_networks([ls.NOT, ""]) == ls.NOT
    assert ls.combine_networks([]) == ls.NOT


def test_balanced_accuracy_always_not_is_a_third():
    pairs = [(ls.UP, ls.NOT)] * 3 + [(ls.DOWN, None)] * 5 + [(ls.NOT, ls.CONFLICT)] * 40
    assert ls.balanced_accuracy(pairs) == pytest.approx(1 / 3)
    perfect = [(m, m) for m, _ in pairs]
    assert ls.balanced_accuracy(perfect) == pytest.approx(1.0)


def test_balanced_accuracy_weights_classes_equally():
    pairs = [(ls.UP, ls.UP)] + [(ls.NOT, ls.NOT)] * 99 + [(ls.DOWN, ls.NOT)]
    assert ls.balanced_accuracy(pairs) == pytest.approx((1 + 1 + 0) / 3)


def test_bootstrap_resamples_proteins_and_reads_the_interval():
    rng = random.Random(1)
    cases = []
    for g in range(40):
        for r in range(10):
            m = rng.choice([ls.UP, ls.DOWN, ls.NOT, ls.NOT])
            cases.append({"gene": f"G{g}", "measured": m, "good": m, "flat": ls.NOT})
    point, lo, hi = ls.bootstrap_difference(cases, "good", "flat", n=500)
    assert point == pytest.approx(2 / 3)
    assert lo > 0 and ls.verdict(lo, hi) == "gain"
    point, lo, hi = ls.bootstrap_difference(cases, "flat", "flat", n=200)
    assert (point, lo, hi) == (0, 0, 0) and ls.verdict(lo, hi) == "unresolved"
    assert ls.verdict(-0.2, -0.01) == "loss"
