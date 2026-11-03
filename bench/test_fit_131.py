"""specs/131: the surrogate fixed point and the implicit gradient of fit_131.py, on a small network with a loop."""
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fit_131 as f  # noqa: E402
import surrogate_131 as sg  # noqa: E402


def _edge(p, c, pos=True, and_=True, t="input"):
    return {"parent_uuid": p, "child_uuid": c, "is_positive": pos, "is_and": and_, "edge_type": t}


NET = {"nodes": [{"uuid": u, "reactome_id": None} for u in "PABCDR"],
       "edges": [_edge("P", "A"), _edge("A", "B"), _edge("B", "C", and_=False, t="output"),
                 _edge("D", "C", and_=False, t="output"), _edge("C", "B", pos=False, and_=False, t="regulator"),
                 _edge("C", "R"), _edge("A", "D")],
       "cofactor_stids": [], "drug_stids": [], "pathogen_stids": []}
OBS = {"P": [80.0, 1.0]}


def test_surrogate_identity_matches_net_forward():
    acts, conv = sg.solve(NET, OBS)
    net = f.Net(NET)
    fixed, s = net.pins(OBS)
    d, conv2 = net.forward(net.matrix(np.zeros(len(net.a))), fixed, s)
    assert conv and conv2
    for u, i in net.idx.items():
        assert abs(acts[u] - 0.01 * 2 ** d[i]) < 1e-12


def test_or_weights_and_signs():
    net = f.Net(NET)
    a = {k[:2]: v for k, v in zip(net.keys, net.a)}
    assert a[("B", "C")] == 0.5 and a[("D", "C")] == 0.5 and a[("C", "B")] == -1.0 and a[("A", "B")] == 1.0


def test_gradient_matches_finite_differences():
    net = f.Net(NET)
    fixed, s = net.pins(OBS)
    rows = [(np.array([net.idx["R"]]), "up"), (np.array([net.idx["B"]]), "down")]
    cw = {"up": 0.5, "down": 0.5, "not": 0.5}
    rng = np.random.default_rng(0)
    theta = rng.normal(0, 0.1, len(net.a))
    s = s * 0.02          # keep every node inside the clip so the derivative exists everywhere
    loss, grad, conv, _ = f.solve_loss_grad(net, theta, fixed, s, rows, cw)
    assert conv
    for k in range(len(theta)):
        e = np.zeros_like(theta)
        e[k] = 1e-6
        lp = f.solve_loss_grad(net, theta + e, fixed, s, rows, cw)[0]
        lm = f.solve_loss_grad(net, theta - e, fixed, s, rows, cw)[0]
        assert abs((lp - lm) / 2e-6 - grad[k]) < 1e-4 * max(1.0, abs(grad[k])), k


def test_route_edges():
    net = f.Net(NET)
    keys = {net.keys[k][:2] for k in net.route_edges([net.idx["P"]], [net.idx["R"]])}
    assert ("P", "A") in keys and ("C", "R") in keys and ("C", "B") in keys


def test_bootstrap_and_balanced_accuracy():
    fitted = [("g1", "up", "up"), ("g1", "not", "not"), ("g2", "down", "down")]
    base = [("g1", "up", "not"), ("g1", "not", "not"), ("g2", "down", "not")]
    point, low, high = f.bootstrap(fitted, base, n=200)
    assert point == 1.0 - 1 / 3 and low <= point <= high
