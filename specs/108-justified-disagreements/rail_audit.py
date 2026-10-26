"""specs/108: which cases are rail-dependent, i.e. the same perturbation at a mild dose
(knockout 0.5x, overexpression 2x) drives the readout to <= 0.01x or >= 99x.

    python3 rail_audit.py DEFAULTS_RESULTS_DIR MILD_RESULTS_DIR [PATHWAY:GENE ...]

Each directory holds `curator_cases.tsv` and `experimental_cases.tsv`. Counts are over cases
valid in both runs; the optional PATHWAY:GENE arguments list those cases one by one.
"""
import csv
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bench" / "analysis"))
from holdout_report import TUNING_PATHWAYS  # noqa: E402

RAIL_LO, RAIL_HI = 0.01, 99.0


def load(path):
    with open(path, newline="") as fh:
        return {(r["pathway"], r["gene"], r["direction"], r["key_output"]): r
                for r in csv.DictReader(fh, delimiter="\t") if r["valid"] == "1" and r["expected"] in "012"}


def railed(r):
    v = float(r["pred_ui"])
    return v <= RAIL_LO or v >= RAIL_HI


def main():
    base, mild, picks = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3:]
    for axis in ("curator", "experimental"):
        a, m = load(base / f"{axis}_cases.tsv"), load(mild / f"{axis}_cases.tsv")
        keys = a.keys() & m.keys()
        print(f"== {axis}: {len(keys)} cases valid in both")
        for split in ("held-out", "tuning", "all"):
            ks = [k for k in keys if split == "all" or (k[0] in TUNING_PATHWAYS) == (split == "tuning")]
            right = [k for k in ks if a[k]["predicted"] == a[k]["expected"]]
            wrong = [k for k in ks if a[k]["predicted"] != a[k]["expected"]]
            rr = [k for k in right if railed(m[k])]
            rw = [k for k in wrong if railed(m[k])]
            print(f"  {split:9s} right {len(right):6d}, rail-dependent {len(rr):5d} ({len(rr) / max(1, len(right)):.2%}); "
                  f"wrong {len(wrong):5d}, rail-dependent {len(rw):5d} ({len(rw) / max(1, len(wrong)):.2%})")
            if split == "all":
                by = Counter(k[0] for k in rr)
                print(f"  rail-dependent right cases in {len(by)} pathways; largest: "
                      + ", ".join(f"{p} {n}" for p, n in by.most_common(8)))
                ho = Counter(k[0] for k in rr if k[0] not in TUNING_PATHWAYS)
                print(f"  held-out only: {len(ho)} pathways")
        for pick in picks:
            p, g = pick.split(":")
            for k in sorted(keys):
                if k[0] == p and k[1] == g:
                    print(f"  {pick} dir {k[2]} readout {k[3]}: defaults {a[k]['pred_ui']} (pred {a[k]['predicted']}, "
                          f"expected {a[k]['expected']}); mild {m[k]['pred_ui']} -> {'RAIL' if railed(m[k]) else 'no rail'}")


if __name__ == "__main__":
    main()
