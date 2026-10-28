import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import train_116 as t  # noqa: E402


def test_bins():
    assert [t.fold_bin(x) for x in ("0", "0.05", "0.3", "0.6", "1", "1.6", "4", "20")] == ["zero", 0, 1, 2, 3, 4, 5, 6]
    assert [t.hop_band(x) for x in ("", "1", "5", "6", "20", "21", "41")] == ["none", 0, 0, 1, 2, 3, 4]


def test_metrics():
    pairs = [("down", "down"), ("not", "down"), ("not", "not"), ("not", "not")]
    assert t.balanced_accuracy(pairs) == 0.75          # down recall 0.5, not recall 1
    assert t.balanced_accuracy([("not", m) for m in ("down", "not", "up")]) == 1 / 3


def test_fit_balanced_rule_and_fallbacks():
    # 30 'not' and 10 'down' in one cell: 'down' wins on share of its class (10/10 > 30/60)
    train = [({"fold": 1, "band": 0}, "not")] * 30 + [({"fold": 1, "band": 0}, "down")] * 10 \
        + [({"fold": 3, "band": 1}, "not")] * 30
    model = t.fit(train, ("fold", "band"))
    assert t.predict(model, {"fold": 1, "band": 0}, ("fold", "band")) == "down"
    assert t.predict(model, {"fold": 3, "band": 1}, ("fold", "band")) == "not"
    # unseen cell in a known band: band fallback; unknown band: NOT
    assert t.predict(model, {"fold": 6, "band": 0}, ("fold", "band")) == "down"
    assert t.predict(model, {"fold": 6, "band": 4}, ("fold", "band")) == "not"


def test_lodo_never_trains_on_the_held_out_drug():
    rows = [{"fold": 1, "band": 0, "measured": "down", "model": "down", "reach": "down", "group": g}
            for g in ("A", "B", "A+B", "C")]
    seen = []
    orig = t.evaluate
    t.evaluate = lambda tr, te: seen.append(({d for r in tr for d in r["drugs"]}, {d for r in te for d in r["drugs"]})) or []
    try:
        t.lodo(rows, "korkut")
    finally:
        t.evaluate = orig
    for trained, tested in seen:
        assert not trained & tested
    assert [len(te) for _, te in seen] == [1, 1, 1]
