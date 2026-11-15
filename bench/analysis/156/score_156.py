#!/usr/bin/env python3
"""specs/156 FR-004: score the registered conditioning arm.

Input: the arm TSV(s) written by arm_156.py (shards concatenate). Per
registered weight, balanced accuracy over {down, not, up} of the conditioned
held-out calls (`predicted`) against the measured truth (`measured`), beside
the unconditioned record (`baseline_predicted`) on EXACTLY the same rows.
The gate (spec.md FR-004, registered before the run): at some weight the
balanced-accuracy delta's 95% bootstrap interval over DRUG SETS (10,000
resamples, seed 156) excludes zero on the positive side.

Pure python, no fitting, nothing chosen from the data.
"""
import collections
import csv
import random
import sys

CLASSES = ("down", "not", "up")
N_BOOT = 10_000
SEED = 156


def balanced(rows, col):
    per = {}
    for c in CLASSES:
        sub = [r for r in rows if r["measured"] == c]
        if sub:
            per[c] = sum(r[col] == c for r in sub) / len(sub)
    return sum(per.values()) / len(per) if per else float("nan"), per


def main():
    rows = []
    for path in sys.argv[1:]:
        rows.extend(csv.DictReader(open(path), delimiter="\t"))
    if not rows:
        raise SystemExit("usage: score_156.py ARM_TSV [ARM_TSV ...]")

    by_w = collections.defaultdict(list)
    for r in rows:
        by_w[r["weight"]].append(r)

    rng = random.Random(SEED)
    verdicts = []
    for w in sorted(by_w, key=float):
        sub = by_w[w]
        sets = sorted({r["drug_set"] for r in sub})
        by_set = collections.defaultdict(list)
        for r in sub:
            by_set[r["drug_set"]].append(r)
        nonconv = sum(r["converged"] == "0" for r in sub)

        arm_ba, arm_per = balanced(sub, "predicted")
        base_ba, base_per = balanced(sub, "baseline_predicted")
        arm_acc = sum(r["predicted"] == r["measured"] for r in sub) / len(sub)
        base_acc = sum(r["baseline_predicted"] == r["measured"] for r in sub) / len(sub)

        deltas = []
        for _ in range(N_BOOT):
            picked = [by_set[rng.choice(sets)] for _ in sets]
            flat = [r for g in picked for r in g]
            a, _ = balanced(flat, "predicted")
            b, _ = balanced(flat, "baseline_predicted")
            deltas.append(a - b)
        deltas.sort()
        lo = deltas[int(0.025 * N_BOOT)]
        hi = deltas[int(0.975 * N_BOOT) - 1]
        delta = arm_ba - base_ba
        passes = lo > 0
        verdicts.append(passes)

        print(f"weight {w}: {len(sub)} held-out rows, {len(sets)} drug sets, "
              f"{nonconv} non-converged")
        print(f"  balanced accuracy: conditioned {arm_ba:.4f} "
              f"vs unconditioned {base_ba:.4f}  delta {delta:+.4f} "
              f"[{lo:+.4f}, {hi:+.4f}] 95% boot over drug sets"
              f"  -> {'PASS' if passes else 'fail'}")
        print(f"  raw accuracy: conditioned {arm_acc:.4f} vs {base_acc:.4f}")
        for c in CLASSES:
            if c in arm_per:
                n_c = sum(r["measured"] == c for r in sub)
                print(f"    recall {c:>4} (n={n_c}): "
                      f"{arm_per[c]:.4f} vs {base_per.get(c, float('nan')):.4f}")

    print(f"\nGATE: {'PASS' if any(verdicts) else 'FAIL'} "
          f"(some registered weight beats the record with the interval "
          f"excluding zero: {any(verdicts)})")


if __name__ == "__main__":
    main()
