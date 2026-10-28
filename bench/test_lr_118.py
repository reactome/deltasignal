import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import lr_118 as t  # noqa: E402


def test_fit_recovers_r_and_c_exactly():
    rng = np.random.default_rng(0)
    D, F = rng.normal(size=(30, 3)), rng.normal(size=(30, 4))
    R = rng.normal(size=(3, 4))
    Rh, c = t.fit(D, D @ R + 0.7 * F, F)
    assert np.allclose(Rh, R) and abs(c - 0.7) < 1e-9
    Rh, c = t.fit(D, D @ R + F * np.array([0.1, 0.2, 0.3, 0.4]), F, per_response=True)
    assert np.allclose(c, [0.1, 0.2, 0.3, 0.4]) and np.allclose(Rh, R)


def test_lodo_never_fits_on_the_held_out_drug():
    sets = [(d,) for d in t.ACTIVITY] + [tuple(sorted((t.ACTIVITY[0], t.ACTIVITY[k]))) for k in range(1, 12)]
    D = np.array([[1.0 if c in s else 0.0 for c in t.ACTIVITY] for s in sets])
    X = np.random.default_rng(1).normal(size=(len(sets), t.N_RESP))
    seen = []
    orig = t.fit
    t.fit = lambda Dtr, Xtr, F=None, per=False: seen.append(Xtr.shape[0]) or orig(Dtr, Xtr, F, per)
    try:
        res, _ = t.lodo(X, D, sets, {})
    finally:
        t.fit = orig
    # drug 0 is in 12 conditions (train 11); every other drug in 2 (train 21)
    assert seen == [11] + [21] * 11
    # LR predicts nothing for a drug given alone: its column is dropped
    assert res["lr"][t.ACTIVITY[5]][1] == 0.0


def test_corr_of_constant_is_zero():
    assert t.corr(np.zeros(5), np.arange(5.0)) == 0.0
