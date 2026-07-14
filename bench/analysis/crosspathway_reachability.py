#!/usr/bin/env python3
"""Quantify the cross-pathway (joint-solving) ceiling.

For every `no_path` failure (readout unreachable from the perturbed gene within
its single pathway), check whether the readout IS reachable once all benchmark
pathways are unioned into one cross-pathway graph. Reachability is computed at
the STABLE_ID level (UUIDs are per-pathway; the same entity has different UUIDs
in different pathways, so cross-pathway connectivity only exists at stid level).

Buckets each no_path case:
  - stid_single : reachable at stid level within its own pathway (UUID-level
                  no_path but stid-level path exists → a within-pathway silo /
                  positional-decomposition artifact, not truly cross-pathway)
  - joint_only  : NOT reachable single, but reachable in the union graph
                  → RECOVERABLE by cross-pathway / joint solving
  - unreachable : not reachable even in the union → genuine Reactome gap

Usage:  python3 bench/analysis/crosspathway_reachability.py
Reads the case dumps /tmp/dsX_{cur,exp}_cases.tsv and /tmp/gene_to_stids.json.
"""
import csv, json, os
from collections import defaultdict, deque
from pathlib import Path

CATALOG = Path(os.environ.get(
    "PATHWAY_CATALOG", str(Path.home() / "gitroot" / "logic-network-generator" / "output")))
EXP_DIRS = [
    "Transcriptional_Regulation_by_TP53_R-HSA-3700989",
    "Mitotic_G1-G1_S_phases_R-HSA-453279", "S_Phase_R-HSA-69242",
    "Signaling_by_ERBB2_R-HSA-1227986", "Cell_Cycle_Checkpoints_R-HSA-69620",
    "Signaling_by_WNT_R-HSA-195721", "PIP3_activates_AKT_signaling_R-HSA-1257604",
    "Mitotic_Prophase_R-HSA-68875",
]
NAME2DIR = {d.split("_R-HSA")[0]: d for d in EXP_DIRS}


def stid_adj(dirname):
    """stid -> set(stid) directed edges for one pathway."""
    d = CATALOG / dirname
    if not d.is_dir():
        return defaultdict(set)
    u2s = {}
    with open(d / "stid_to_uuid_mapping.csv") as f:
        for r in csv.DictReader(f):
            u2s[r["uuid"]] = r["stable_id"]
    adj = defaultdict(set)
    with open(d / "logic_network.csv") as f:
        for r in csv.DictReader(f):
            s, t = u2s.get(r["source_id"]), u2s.get(r["target_id"])
            if s and t:
                adj[s].add(t)
    return adj


def reachable(adj, sources, target, cap=200000):
    seen = set(sources)
    q = deque(sources)
    steps = 0
    while q and steps < cap:
        u = q.popleft(); steps += 1
        if u == target:
            return True
        for v in adj.get(u, ()):
            if v not in seen:
                seen.add(v); q.append(v)
    return target in seen


def main():
    gene2stids = json.load(open("/tmp/gene_to_stids.json"))
    per_pw = {name: stid_adj(d) for name, d in NAME2DIR.items()}
    # combined union graph
    combined = defaultdict(set)
    for adj in per_pw.values():
        for s, ts in adj.items():
            combined[s] |= ts
    print(f"combined cross-pathway graph: {len(combined)} stid nodes\n")

    for dump, lbl in [("/tmp/dsX_cur_cases.tsv", "CURATOR"),
                      ("/tmp/dsX_exp_cases.tsv", "EXPERIMENTAL")]:
        if not os.path.exists(dump):
            continue
        rows = [r for r in csv.DictReader(open(dump), delimiter="\t")
                if r["category"] == "no_path"]
        buckets = defaultdict(lambda: defaultdict(int))
        for r in rows:
            pw = r["pathway"]
            gstids = gene2stids.get(r["gene"], [])
            # readout stid: try R-HSA- then R-ALL-
            ko = r["key_output"]
            cands = [f"R-HSA-{ko}", f"R-ALL-{ko}", ko]
            single = per_pw.get(pw, defaultdict(set))
            # target present in combined?
            tgt = next((c for c in cands if c in combined or c in single), cands[0])
            r_single = reachable(single, [g for g in gstids if g in single], tgt)
            r_joint = reachable(combined, [g for g in gstids if g in combined], tgt)
            if r_single:
                buckets[pw]["stid_single(silo)"] += 1
            elif r_joint:
                buckets[pw]["joint_only(RECOVERABLE)"] += 1
            else:
                buckets[pw]["unreachable(gap)"] += 1
        print(f"########## {lbl}  no_path cases ##########")
        tot = defaultdict(int)
        for pw in sorted(buckets):
            b = buckets[pw]; n = sum(b.values())
            for k, v in b.items(): tot[k] += v
            parts = "  ".join(f"{k}={v}" for k, v in sorted(b.items()))
            print(f"  {pw[:32]:32} n={n:4}  {parts}")
        print(f"  {'TOTAL':32}        " + "  ".join(f"{k}={v}" for k, v in sorted(tot.items())))
        print()


if __name__ == "__main__":
    main()
