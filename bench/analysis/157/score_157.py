#!/usr/bin/env python3
"""specs/157: label robust/fragile per the registered R2 rule and run the gate.

Input: the ensemble runs' calls_*.tsv (calls_canonical.tsv required), the
registry, and the 98 record's case TSVs for truth and the canonical
self-check. Labels (R2, registered before any run): per (axis, case, readout),
a variant CONTRIBUTES only when its solve converged; ROBUST = every
contributing call equals the canonical run's call; FRAGILE = any differs; a
pair with fewer than 80% of variants contributing is EXCLUDED and counted.
The gate (FR-003): on at least one axis, accuracy(robust) > accuracy(fragile)
with a 95% bootstrap interval over the axis's independent units (perturbations
on curator/experimental, drug sets on korkut), 10k resamples, seed 157,
excluding zero. Everything reported whatever it says; RAF read beside the
curator headline per the standing rule.
"""
import collections
import csv
import random
import sys
from pathlib import Path

N_BOOT = 10_000
SEED = 157
CONTRIB_MIN = 0.80
CODE = {"0": "down", "1": "not", "2": "up"}
RAF = "RAF_MAP_kinase_cascade"


def load_calls(path):
    calls = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        key = (r["axis"], r["pathway"], r["case_a"], r["case_b"], r["readout"])
        calls[key] = (r["call"], r["converged"] == "1")
    return calls


def truth_maps(curator_tsv, experimental_tsv, korkut_tsv):
    truth, canon_rec = {}, {}
    for axis, path in (("curator", curator_tsv), ("experimental", experimental_tsv)):
        for r in csv.DictReader(open(path), delimiter="\t"):
            if r.get("valid") != "1":
                continue
            key = (axis, r["pathway"], r["gene"], r["direction"], r["key_output"])
            truth[key] = CODE[r["expected"]]
            canon_rec[key] = CODE[r["predicted"]]
    for r in csv.DictReader(open(korkut_tsv), delimiter="\t"):
        if r["status"] != "scored":
            continue
        key = ("korkut", r["pathway_id"], r["drug_set"], "", r["readout"])
        truth[key] = r["measured"]
        canon_rec[key] = r["predicted"]
    return truth, canon_rec


def unit_of(key):
    axis = key[0]
    return (key[1], key[2], key[3]) if axis in ("curator", "experimental") else key[2]


def gate(pairs, rng):
    """pairs: list of (unit, robust: bool, correct: bool) -> stats dict."""
    by_unit = collections.defaultdict(list)
    for unit, robust, correct in pairs:
        by_unit[unit].append((robust, correct))
    units = sorted(by_unit)

    def acc(rows, want_robust):
        sub = [c for (r, c) in rows if r == want_robust]
        return (sum(sub) / len(sub), len(sub)) if sub else (float("nan"), 0)

    flat = [rc for u in units for rc in by_unit[u]]
    rob, n_rob = acc(flat, True)
    fra, n_fra = acc(flat, False)
    deltas = []
    for _ in range(N_BOOT):
        picked = [by_unit[rng.choice(units)] for _ in units]
        rows = [rc for g in picked for rc in g]
        a, na = acc(rows, True)
        b, nb = acc(rows, False)
        deltas.append((a - b) if na and nb else 0.0)
    deltas.sort()
    lo = deltas[int(0.025 * N_BOOT)]
    hi = deltas[int(0.975 * N_BOOT) - 1]
    return dict(robust_acc=rob, n_robust=n_rob, fragile_acc=fra, n_fragile=n_fra,
                delta=rob - fra, lo=lo, hi=hi, passes=lo > 0, units=len(units))


def main():
    out_dir = Path(sys.argv[1])
    curator_tsv, experimental_tsv, korkut_tsv = sys.argv[2], sys.argv[3], sys.argv[4]
    canonical = load_calls(out_dir / "calls_canonical.tsv")
    variants = {}
    for p in sorted(out_dir.glob("calls_*.tsv")):
        name = p.stem[len("calls_"):]
        if name != "canonical":
            variants[name] = load_calls(p)
    if not variants:
        raise SystemExit("no variant call files found")
    truth, canon_rec = truth_maps(curator_tsv, experimental_tsv, korkut_tsv)

    # self-check: the canonical run's calls vs the recorded predicted
    both = [k for k in canonical if k in canon_rec and canonical[k][0]]
    mism = [k for k in both if canonical[k][0] != canon_rec[k]]
    print(f"self-check: canonical run reproduces the recorded call on "
          f"{len(both) - len(mism)} of {len(both)} truth-joined pairs "
          f"({len(mism)} mismatches)")

    labels = {}
    excluded = 0
    flips_by_variant = collections.Counter()
    for key, (ccall, _) in canonical.items():
        if key not in truth or not ccall:
            continue
        contributing = 0
        fragile = False
        for name, calls in variants.items():
            vcall, conv = calls.get(key, ("", False))
            if not conv or not vcall:
                continue
            contributing += 1
            if vcall != ccall:
                fragile = True
                flips_by_variant[name] += 1
        if contributing < CONTRIB_MIN * len(variants):
            excluded += 1
            continue
        labels[key] = fragile
    print(f"{len(labels)} labelled pairs ({excluded} excluded below the "
          f"{CONTRIB_MIN:.0%} contribution floor) over {len(variants)} variants")

    rng = random.Random(SEED)
    any_pass = False
    for axis in ("curator", "experimental", "korkut"):
        pairs = [(unit_of(k), not fragile, canonical[k][0] == truth[k])
                 for k, fragile in labels.items() if k[0] == axis]
        if not pairs:
            continue
        g = gate(pairs, rng)
        any_pass = any_pass or g["passes"]
        frag_share = g["n_fragile"] / (g["n_robust"] + g["n_fragile"])
        print(f"{axis}: robust {g['robust_acc']:.4f} (n={g['n_robust']}) vs "
              f"fragile {g['fragile_acc']:.4f} (n={g['n_fragile']}, "
              f"{frag_share:.1%} of pairs) delta {g['delta']:+.4f} "
              f"[{g['lo']:+.4f}, {g['hi']:+.4f}] over {g['units']} units "
              f"-> {'PASS' if g['passes'] else 'fail'}")
        if axis == "curator":
            nr = [(u, r, c) for (u, r, c) in pairs if u[0] != RAF]
            if len(nr) < len(pairs):
                g2 = gate(nr, rng)
                print(f"  (RAF set aside: robust {g2['robust_acc']:.4f} vs "
                      f"fragile {g2['fragile_acc']:.4f} delta {g2['delta']:+.4f} "
                      f"[{g2['lo']:+.4f}, {g2['hi']:+.4f}])")

    print("\nflips by variant (a pair may flip under several):")
    for name, n in flips_by_variant.most_common():
        print(f"  {name}: {n}")
    print(f"\nGATE: {'PASS' if any_pass else 'FAIL'}")

    with open(out_dir / "labels.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["axis", "pathway", "case_a", "case_b", "readout", "fragile",
                    "canonical_call", "truth"])
        for k, fragile in sorted(labels.items()):
            w.writerow([*k, int(fragile), canonical[k][0], truth[k]])


if __name__ == "__main__":
    main()
