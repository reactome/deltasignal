"""specs/107: score behaviour-check runs (B2 convergence, B3 rails, B4 dose
monotonicity) per arm, on the pairs every arm solved.

    python3 score_behaviour.py CONTROL.tsv ARM.tsv [ARM.tsv ...]

Each file is `behaviour.jl run` output. A solve is (pathway, root, dose); a pair is
(pathway, root, readout). Rails: a readout >= 99x or <= 0.01x under a 0.5x or 2x
root. Monotone: the readout's folds over root doses 0, 0.5, 1 (the unperturbed
solve), 2, 10 never rise and fall (relative tolerance 1e-6).
"""
import csv
import sys
from collections import defaultdict

TOL = 1e-6


def load(path):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    solves, vals, base = {}, {}, {}
    for r in rows:
        dose = float(r["dose"])
        if r["root"] == "":
            base[(r["pathway"], r["readout"])] = float(r["fold"])
            solves[(r["pathway"], "", 1.0)] = r["converged"] == "true"
            continue
        solves[(r["pathway"], r["root"], dose)] = r["converged"] == "true"
        vals[(r["pathway"], r["root"], r["readout"], dose)] = float(r["fold"])
    return solves, vals, base


def monotone(xs):
    up = all(b >= a - TOL * max(1.0, abs(a)) for a, b in zip(xs, xs[1:]))
    down = all(b <= a + TOL * max(1.0, abs(a)) for a, b in zip(xs, xs[1:]))
    return up or down


def score(path, keys=None):
    solves, vals, base = load(path)
    pairs = sorted({(p, r, o) for (p, r, o, _) in vals})
    if keys is not None:
        pairs = [k for k in pairs if k in keys]
    conv = [v for k, v in solves.items()]
    rails = total = mono = mono_n = 0
    by_pw = defaultdict(lambda: [0, 0, 0, 0])
    for (p, r, o) in pairs:
        for d in (0.5, 2.0):
            x = vals.get((p, r, o, d))
            if x is None:
                continue
            total += 1
            if x >= 99.0 or x <= 0.01:
                rails += 1
                by_pw[p][0] += 1
            by_pw[p][1] += 1
        seq = [vals.get((p, r, o, 0.0)), vals.get((p, r, o, 0.5)), base.get((p, o)),
               vals.get((p, r, o, 2.0)), vals.get((p, r, o, 10.0))]
        if None in seq:
            continue
        mono_n += 1
        ok = monotone(seq)
        mono += ok
        by_pw[p][2] += ok
        by_pw[p][3] += 1
    return {"solves": len(conv), "converged": sum(conv), "pairs": len(pairs),
            "rails": rails, "rail_checks": total, "monotone": mono, "monotone_n": mono_n}, by_pw, set(pairs)


def main():
    files = sys.argv[1:]
    keys = None
    for f in files:
        _, _, ks = score(f)
        keys = ks if keys is None else keys & ks
    ctl = None
    for f in files:
        s, by_pw, _ = score(f, keys)
        line = (f"{f}: solves {s['solves']}, converged {s['converged']} ({s['converged']/s['solves']:.1%}); "
                f"rails {s['rails']}/{s['rail_checks']} ({s['rails']/max(1, s['rail_checks']):.2%}); "
                f"monotone {s['monotone']}/{s['monotone_n']} ({s['monotone']/max(1, s['monotone_n']):.2%})")
        print(line)
        if ctl is None:
            ctl = by_pw
            continue
        for p in sorted(set(ctl) | set(by_pw)):
            c, a = ctl.get(p, [0, 0, 0, 0]), by_pw.get(p, [0, 0, 0, 0])
            if c[0] != a[0] or c[2] != a[2]:
                print(f"    {p}: rails {c[0]} -> {a[0]} of {a[1]}; monotone {c[2]} -> {a[2]} of {a[3]}")


if __name__ == "__main__":
    main()
