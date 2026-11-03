#!/usr/bin/env python3
"""specs/131 S2: fit a strength per edge of the surrogate on korkut, cross-validated by drug set.

Each edge's entry in A is scaled by W = exp(theta). A row's value is z = log2(mean over its readout nodes of 2^delta);
its class probabilities are softmax(k (z - u), 0, k (d - z)) for up / not / down, u = log2 1.15, d = log2 0.85, and the
loss is cross-entropy weighted so each measured class counts equally (balanced accuracy's weighting), plus
lam * mean(theta^2) over the free edges. Gradients by implicit differentiation at the fixed point over free nodes that
are neither pinned, held nor clipped; solves that do not converge at the current weights are left out of the loss.

Usage: fit_131.py PROBLEMS.pickle OUT.json [--quick]
"""
import json
import math
import multiprocessing as mp
import pickle
import random
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as spla

sys.path.insert(0, str(Path(__file__).resolve().parent))
import surrogate_131 as sg  # noqa: E402

UP, DOWN = math.log2(1.15), math.log2(0.85)
K = 10.0
CLASSES = ("up", "not", "down")
SEED = 131
FOLDS = 5
INNER = 3
LAMS = (0.01, 0.1, 1.0, 10.0)
EPOCHS = 40
LR = 0.05
N_BOOT = 10000


# ---------------------------------------------------------------- one network, prepared once

class Net:
    def __init__(self, network):
        self.idx, A, edges = sg.build(network)
        self.n = len(self.idx)
        self.rows = np.array([e[0] for e in edges], dtype=np.int64)
        self.cols = np.array([e[1] for e in edges], dtype=np.int64)
        self.a = np.array([e[2] for e in edges])
        self.keys = [e[3] for e in edges]
        self.held = {self.idx[u] for u in sg.held_uuids(network) if u in self.idx}
        self.children = defaultdict(list)
        self.parents = defaultdict(list)
        for k, (c, p) in enumerate(zip(self.rows, self.cols)):
            self.children[int(p)].append(int(c))
            self.parents[int(c)].append(int(p))

    def matrix(self, theta):
        return sparse.csr_matrix((self.a * np.exp(theta), (self.rows, self.cols)), shape=(self.n, self.n))

    def pins(self, obs):
        fixed = np.zeros(self.n, dtype=bool)
        s = np.zeros(self.n)
        for u, (ui, conf) in obs.items():
            if u in self.idx and conf > 1e-6:
                i = self.idx[u]
                fixed[i] = True
                s[i] = -sg.CLIP if ui <= 0 else max(-sg.CLIP, min(sg.CLIP, math.log2(ui)))
        for i in self.held:
            fixed[i] = True
        return fixed, s

    def forward(self, M, fixed, s, start=None, sweeps=sg.SWEEPS):
        d = s.copy() if start is None else start.copy()
        d[fixed] = s[fixed]
        for _ in range(sweeps):
            new = np.clip(M @ d, -sg.CLIP, sg.CLIP)
            new[fixed] = s[fixed]
            if np.max(np.abs(new - d)) < sg.TOL:
                return new, True
            d = new
        return d, False

    def route_edges(self, pin_nodes, readout_nodes):
        """Edges with a parent reachable from the pins and a child that reaches a readout node."""
        down, stack = set(pin_nodes), list(pin_nodes)
        while stack:
            for c in self.children[stack.pop()]:
                if c not in down:
                    down.add(c), stack.append(c)
        up, stack = set(readout_nodes), list(readout_nodes)
        while stack:
            for p in self.parents[stack.pop()]:
                if p not in up:
                    up.add(p), stack.append(p)
        return {k for k, (c, p) in enumerate(zip(self.rows, self.cols)) if int(p) in down and int(c) in up}


# ---------------------------------------------------------------- loss and gradient for one solve

def row_value(d, nodes):
    v = 2.0 ** d[nodes]
    return math.log2(v.mean()), v / v.sum()


def probs(z):
    l = np.array([K * (z - UP), 0.0, K * (DOWN - z)])
    l -= l.max()
    e = np.exp(l)
    return e / e.sum()


def call(z):
    return "up" if z >= UP else "down" if z <= DOWN else "not"


def solve_loss_grad(net, theta, fixed, s, rows, cw, start=None):
    """(loss, grad over the network's edges, converged, final delta) for one solve; rows [(nodes, measured)]."""
    M = net.matrix(theta)
    d, conv = net.forward(M, fixed, s, start)
    if not conv:
        return 0.0, np.zeros(len(net.a)), False, d
    loss, g = 0.0, np.zeros(net.n)
    for nodes, m in rows:
        z, w = row_value(d, nodes)
        p = probs(z)
        j = CLASSES.index(m)
        loss -= cw[m] * math.log(max(p[j], 1e-12))
        # dCE/dz for softmax over logits (k(z-u), 0, k(d-z)): sum_c (p_c - y_c) dl_c/dz
        y = np.zeros(3)
        y[j] = 1.0
        dz = cw[m] * ((p[0] - y[0]) * K - (p[2] - y[2]) * K)
        g[nodes] += dz * w
    free = (~fixed) & (np.abs(d) < sg.CLIP - 1e-9)
    U = np.flatnonzero(free)
    grad = np.zeros(len(net.a))
    if len(U) == 0 or not np.any(g[U]):
        return loss, grad, True, d
    pos = -np.ones(net.n, dtype=np.int64)
    pos[U] = np.arange(len(U))
    Muu = M[U][:, U]
    lhs = (sparse.identity(len(U), format="csr") - Muu).T.tocsc()
    try:
        lam = spla.spsolve(lhs, g[U])
    except Exception:
        return loss, grad, True, d
    if not np.all(np.isfinite(lam)):
        return loss, grad, True, d
    ci = pos[net.rows]
    ok = ci >= 0
    grad[ok] = lam[ci[ok]] * net.a[ok] * np.exp(theta[ok]) * d[net.cols[ok]]
    return loss, grad, True, d


# ---------------------------------------------------------------- data and parallel evaluation

_DATA = {}


def _init(problems_path):
    with open(problems_path, "rb") as f:
        prob = pickle.load(f)
    nets = {k: Net(v) for k, v in prob["networks"].items()}
    solves = []
    for sv in prob["solves"]:
        net = nets[sv["pathway"]]
        fixed, s = net.pins(sv["obs"])
        rows = [(np.array([net.idx[u] for u in free if u in net.idx]), m) for _, free, m in sv["rows"]]
        rows = [(nodes, m) for nodes, m in rows if len(nodes) and m in CLASSES]   # a 'conflict' call is not scored
        solves.append({"pathway": sv["pathway"], "drug_set": sv["drug_set"], "fixed": fixed, "s": s,
                       "rows": rows, "names": [n for n, f, _ in sv["rows"]]})
    _DATA.update(nets=nets, solves=solves)


def _work(args):
    i, theta, cw, mode = args
    sv = _DATA["solves"][i]
    net = _DATA["nets"][sv["pathway"]]
    if mode == "grad":
        loss, grad, conv, _ = solve_loss_grad(net, theta, sv["fixed"], sv["s"], sv["rows"], cw)
        return i, loss, grad, conv
    d, conv = net.forward(net.matrix(theta), sv["fixed"], sv["s"])
    return i, [call(row_value(d, nodes)[0]) for nodes, _ in sv["rows"]], None, conv


class Model:
    """theta for every network edge, stored per network."""

    def __init__(self, pool):
        self.pool = pool
        self.theta = {p: np.zeros(len(n.a)) for p, n in _DATA["nets"].items()}

    def calls(self, ids):
        out = {}
        for i, cs, _, _ in self.pool.imap_unordered(_work, [(i, self.theta[_DATA["solves"][i]["pathway"]], None,
                                                              "eval") for i in ids], chunksize=4):
            out[i] = cs
        return out

    def fit(self, ids, lam, free, epochs=EPOCHS):
        counts = defaultdict(int)
        for i in ids:
            for _, m in _DATA["solves"][i]["rows"]:
                counts[m] += 1
        n_rows = sum(counts.values())
        cw = {c: n_rows / (len(counts) * counts[c]) / n_rows for c in counts}   # mean balanced CE
        mask = {p: np.zeros(len(n.a), dtype=bool) for p, n in _DATA["nets"].items()}
        for p, ks in free.items():
            mask[p][list(ks)] = True
        n_free = max(1, sum(int(m.sum()) for m in mask.values()))
        m1 = {p: np.zeros_like(t) for p, t in self.theta.items()}
        m2 = {p: np.zeros_like(t) for p, t in self.theta.items()}
        history = []
        for ep in range(1, epochs + 1):
            grads = {p: np.zeros_like(t) for p, t in self.theta.items()}
            loss, conv = 0.0, 0
            jobs = [(i, self.theta[_DATA["solves"][i]["pathway"]], cw, "grad") for i in ids]
            for i, l, g, c in self.pool.imap_unordered(_work, jobs, chunksize=4):
                p = _DATA["solves"][i]["pathway"]
                loss += l
                grads[p] += g
                conv += c
            reg = 0.0
            for p in self.theta:
                t = self.theta[p]
                reg += float(np.sum(t[mask[p]] ** 2))
                g = np.where(mask[p], grads[p] + 2 * lam * t / n_free, 0.0)
                m1[p] = 0.9 * m1[p] + 0.1 * g
                m2[p] = 0.999 * m2[p] + 0.001 * g * g
                step = LR * (m1[p] / (1 - 0.9 ** ep)) / (np.sqrt(m2[p] / (1 - 0.999 ** ep)) + 1e-8)
                self.theta[p] = np.where(mask[p], t - step, 0.0)
            history.append({"epoch": ep, "loss": loss, "reg": lam * reg / n_free, "converged": conv})
        return history


# ---------------------------------------------------------------- folds, choice, scoring

def assign_folds(groups, k, seed):
    order = sorted(set(groups))
    random.Random(seed).shuffle(order)
    return {g: i % k for i, g in enumerate(order)}


def balanced_accuracy(pairs):
    hit, tot = defaultdict(int), defaultdict(int)
    for m, p in pairs:
        tot[m] += 1
        hit[m] += int(m == p)
    present = [c for c in CLASSES if tot[c]]
    return sum(hit[c] / tot[c] for c in present) / len(present) if present else float("nan")


def free_edges(ids):
    out = defaultdict(set)
    for i in ids:
        sv = _DATA["solves"][i]
        net = _DATA["nets"][sv["pathway"]]
        pins = [j for j in np.flatnonzero(sv["fixed"]) if j not in net.held]
        ro = sorted({int(x) for nodes, _ in sv["rows"] for x in nodes})
        out[sv["pathway"]] |= net.route_edges(pins, ro)
    return out


def score(pool, ids, theta=None):
    m = Model(pool)
    if theta is not None:
        m.theta = theta
    calls = m.calls(ids)
    return [(_DATA["solves"][i]["drug_set"], meas, c) for i in ids
            for (_, meas), c in zip(_DATA["solves"][i]["rows"], calls[i])]


def bootstrap(fitted, base, n=N_BOOT, seed=SEED):
    """Balanced accuracy of fitted − base and a 95% interval from resampling drug sets."""
    by = defaultdict(list)
    for (g, m, p), (_, _, q) in zip(fitted, base):
        by[g].append((m, p, q))
    groups = sorted(by)

    def diff(sample):
        return (balanced_accuracy([(m, p) for m, p, _ in sample])
                - balanced_accuracy([(m, q) for m, _, q in sample]))

    rng = random.Random(seed)
    draws = sorted(diff([x for g in (rng.choice(groups) for _ in groups) for x in by[g]]) for _ in range(n))
    return diff([x for g in groups for x in by[g]]), draws[int(0.025 * n)], draws[int(0.975 * n) - 1]


def main():
    problems, out = sys.argv[1:3]
    quick = "--quick" in sys.argv
    epochs = 3 if quick else EPOCHS
    lams = LAMS[:1] if quick else LAMS
    _init(problems)
    solves = _DATA["solves"]
    folds = assign_folds([s["drug_set"] for s in solves], FOLDS, SEED)
    ids_all = list(range(len(solves)))
    report = {"solves": len(solves), "rows": sum(len(s["rows"]) for s in solves), "folds": []}
    cv_scores = []
    with mp.get_context("fork").Pool(8, initializer=_init, initargs=(problems,)) as pool:
        base = score(pool, ids_all)
        for f in range(FOLDS):
            train = [i for i in ids_all if folds[solves[i]["drug_set"]] != f]
            test = [i for i in ids_all if folds[solves[i]["drug_set"]] == f]
            # inner CV for lam over the training drug sets
            inner = assign_folds([solves[i]["drug_set"] for i in train], INNER, SEED + f)
            lam_score = {}
            for lam in lams:
                pairs = []
                for g in range(INNER):
                    itr = [i for i in train if inner[solves[i]["drug_set"]] != g]
                    ite = [i for i in train if inner[solves[i]["drug_set"]] == g]
                    m = Model(pool)
                    m.fit(itr, lam, free_edges(itr), epochs)
                    pairs += [(meas, c) for _, meas, c in score(pool, ite, m.theta)]
                lam_score[lam] = balanced_accuracy(pairs)
            lam = max(lams, key=lambda x: (round(lam_score[x], 12), x))     # ties: the stronger pull to 1
            m = Model(pool)
            hist = m.fit(train, lam, free_edges(train), epochs)
            fold_scores = score(pool, test, m.theta)
            cv_scores += [(i, x) for i, x in zip([i for i in test for _ in solves[i]["rows"]], fold_scores)]
            moved = [float(t) for p in m.theta for t in m.theta[p] if t != 0.0]
            report["folds"].append({"fold": f, "lam_scores": lam_score, "lam": lam, "history": hist[-1],
                                    "free_edges": sum(len(v) for v in free_edges(train).values()),
                                    "theta_moved": len(moved),
                                    "theta_quantiles": np.quantile(moved, [0.05, 0.5, 0.95]).tolist() if moved else []})
            print(json.dumps(report["folds"][-1]), flush=True)
    order = {i: k for k, i in enumerate(ids_all)}
    cv_scores.sort(key=lambda x: order[x[0]])
    fitted = [x for _, x in cv_scores]
    gain = bootstrap(fitted, base)
    report.update(base_balanced_accuracy=balanced_accuracy([(m, c) for _, m, c in base]),
                  cv_balanced_accuracy=balanced_accuracy([(m, c) for _, m, c in fitted]),
                  cv_gain=dict(zip(("point", "low", "high"), gain)), passes=gain[1] > 0)
    with open(out, "w") as fh:
        json.dump(report, fh, indent=1)
    print(json.dumps({k: v for k, v in report.items() if k != "folds"}, indent=1))


if __name__ == "__main__":
    main()
