#!/usr/bin/env python3
"""specs/153 S2: per-arrow lead rates on decryptM, F discipline (specs/148).

Join control and composite-arm decryptm case tables on the case key; on
control-untouched cases (control scored, pred_ui == 1) take the arm's answer;
a bridged change call is attributed to the (centre, drug) arrow from the S1
coverage table; rate = share matching the measured regulation, Clopper-Pearson
95% interval.

Usage: s2_rate.py CONTROL_TSV ARM_TSV S1_COVERAGE_TSV
"""
import csv, sys
from collections import defaultdict
import math

def _beta_ppf(q, a, b):
    """Clopper-Pearson bound by bisection on the regularized incomplete beta."""
    def betainc(a, b, x):
        # continued fraction (Lentz), as in Numerical Recipes betacf
        if x in (0.0, 1.0):
            return x
        lbeta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
        front = math.exp(lbeta + a * math.log(x) + b * math.log(1 - x))
        if x < (a + 1) / (a + b + 2):
            return front * _betacf(a, b, x) / a
        return 1 - front * _betacf(b, a, 1 - x) / b
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if betainc(a, b, mid) < q:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2

def _betacf(a, b, x):
    MAXIT, EPS, FPMIN = 200, 3e-12, 1e-300
    qab, qap, qam = a + b, a + 1, a - 1
    c, d = 1.0, 1 - qab * x / qap
    if abs(d) < FPMIN: d = FPMIN
    d = 1 / d
    h = d
    for m in range(1, MAXIT + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1 + aa * d
        if abs(d) < FPMIN: d = FPMIN
        c = 1 + aa / c
        if abs(c) < FPMIN: c = FPMIN
        d = 1 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1 + aa * d
        if abs(d) < FPMIN: d = FPMIN
        c = 1 + aa / c
        if abs(c) < FPMIN: c = FPMIN
        d = 1 / d
        de = d * c
        h *= de
        if abs(de - 1) < EPS:
            break
    return h

class beta:
    @staticmethod
    def ppf(q, a, b):
        return _beta_ppf(q, a, b)

ctl_p, arm_p, cov_p = sys.argv[1:4]

def load(p):
    return {(r["pathway_id"], r["drug"], r["targets"], r["uniprot"], r["position"]): r
            for r in csv.DictReader(open(p), delimiter="\t")}

ctl, arm = load(ctl_p), load(arm_p)
cov = defaultdict(set)   # (centre, drug) -> arrow label
for r in csv.DictReader(open(cov_p), delimiter="\t"):
    if r["target_reaches_entity"] == "True":
        cov[(r["centre"], r["drug"])].add(f"{r['source_name']} -> {r['entity_name']} -> {r['centre_name']}")

keys = set(ctl) & set(arm)
print(f"cases joined: {len(keys)} (control {len(ctl)}, arm {len(arm)})")

stats = defaultdict(lambda: {"bridged": 0, "right": 0, "false": 0, "inverted": 0})
unattr = 0
for k in sorted(keys):
    c, a = ctl[k], arm[k]
    if c["status"] != "scored" or a["status"] != "scored":
        continue
    if float(c["pred_ui"]) != 1.0:
        continue                       # centre routes on its own: its answer stands
    if a["predicted"] == "not":
        continue                       # no bridged change call
    arrows = cov.get((k[0], k[1]))
    if not arrows:
        unattr += 1
        continue
    if len(arrows) > 1:
        label = " / ".join(sorted(arrows)) + " [AMBIGUOUS]"
    else:
        label = next(iter(arrows))
    s = stats[label]
    s["bridged"] += 1
    m = c["measured"]
    if m == a["predicted"]:
        s["right"] += 1
    elif m == "not":
        s["false"] += 1
    else:
        s["inverted"] += 1

print(f"bridged change calls with no covered arrow for (centre, drug): {unattr}")
for label, s in sorted(stats.items()):
    n, r = s["bridged"], s["right"]
    lo = beta.ppf(0.025, r, n - r + 1) if r > 0 else 0.0
    hi = beta.ppf(0.975, r + 1, n - r) if r < n else 1.0
    print(f"\n{label}")
    print(f"  bridged change calls {n}: right {r}, false {s['false']}, inverted {s['inverted']}")
    print(f"  precision {r/n:.3f} [{lo:.3f}, {hi:.3f}]" + ("  (n < 20: reported, stays unrated)" if n < 20 else ""))
