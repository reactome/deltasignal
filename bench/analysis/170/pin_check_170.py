#!/usr/bin/env python3
"""specs/170 FR-003: for every curated node of every pathway in a literature arm, is it a root under
DS_PIN_CURATED_ROOT exactly when the curated network makes it one? Every difference is listed with its reason.
Also counts the curated roots a literature edge unroots without the rule (the specs/168 mechanism).

Usage: pin_check_170.py ARM_DIR CURATED_DIR OUT_TSV      (exit 1 on an unexplained difference)
"""
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
import benchmark_vs_mpbiopath as bm  # noqa: E402
from network_files import open_network  # noqa: E402


def check(arm, curated):
    out, unrooted = [], 0
    for pdir in sorted(p for p in arm.iterdir() if (p / "logic_network.csv").exists()):
        cdir = curated / pdir.name
        cnodes = bm._node_uuids(cdir)
        with open_network(cdir) as f:
            crows = list(csv.DictReader(f))
        with open_network(pdir) as f:
            srows = list(csv.DictReader(f))
        c_in = bm.root_indegree(crows)
        bm.PIN_CURATED_ROOT = str(curated)
        try:
            f_in = bm.root_indegree(bm.curated_pin_rows(pdir, srows))
        finally:
            bm.PIN_CURATED_ROOT = ""
        s_in = bm.root_indegree(srows)
        c_out = {r["source_id"] for r in crows if r.get("edge_type") != bm.RECYCLE_EDGE}
        f_out = {r["source_id"] for r in srows if r.get("edge_type") != bm.RECYCLE_EDGE}
        for u in sorted(cnodes):
            c_root, f_root = c_in.get(u, 0) == 0, f_in.get(u, 0) == 0
            unrooted += c_root and s_in.get(u, 0) > 0
            if c_root == f_root:
                continue
            why = ("gains a literature out-edge (a depletion-fed dead end in the curated network)"
                   if u not in c_out and u in f_out else "UNEXPLAINED")
            out.append({"pathway": pdir.name, "uuid": u, "curated_root": int(c_root), "rule_root": int(f_root),
                        "reason": why})
    return out, unrooted


def main():
    arm, curated, out_tsv = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    diffs, unrooted = check(arm, curated)
    with open(out_tsv, "w", newline="") as f:
        w = csv.DictWriter(f, ["pathway", "uuid", "curated_root", "rule_root", "reason"], delimiter="\t")
        w.writeheader()
        w.writerows(diffs)
    bad = [d for d in diffs if d["reason"] == "UNEXPLAINED"]
    print(f"{arm.name}: curated roots a literature edge unroots without the rule: {unrooted}; "
          f"root differences under the rule: {len(diffs)} ({len(bad)} unexplained)")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
