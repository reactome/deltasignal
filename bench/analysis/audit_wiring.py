#!/usr/bin/env python3
"""Audit generator wiring hazards #7 (boundary-def mismatch) and #8 (arbitrary
regulator/catalyst UUID reuse) on the generated logic networks.

CSV-only (no Neo4j) — safe to run alongside the benchmark. Quantifies whether
the latent hazards actually manifest at material scale before deciding to touch
the generator.

#8: regulator/catalyst edges reuse the FIRST-inserted UUID for a stId
    (stid_to_existing_uuid). If a stId has several position-aware UUIDs, that
    choice can attach the regulator to a non-functional occurrence. We measure,
    per catalyst/regulator edge whose source stId has >=2 UUIDs, whether the
    chosen source UUID is the HIGHEST-degree occurrence of that stId (sensible)
    or a lower/zero-degree one (arbitrary/suspect).

#7: Phase-1 boundary UUID sharing uses a GLOBAL stId set-difference
    (root = produced-by-nobody). An entity that is a root input at one position
    but produced elsewhere is excluded. We count stIds that are BOTH produced
    (some UUID is an edge target) and a root (some UUID is never a target).

Usage:
  python3 bench/analysis/audit_wiring.py            # 9 exp pathways
  python3 bench/analysis/audit_wiring.py TP53 WNT
"""
import csv, os, sys
from collections import defaultdict
from pathlib import Path

CATALOG = Path(os.environ.get(
    "PATHWAY_CATALOG", str(Path.home() / "gitroot" / "logic-network-generator" / "output")))

EXP_DIRS = [
    "Transcriptional_Regulation_by_TP53_R-HSA-3700989",
    "Mitotic_G1-G1_S_phases_R-HSA-453279",
    "S_Phase_R-HSA-69242",
    "Signaling_by_ERBB2_R-HSA-1227986",
    "Cell_Cycle_Checkpoints_R-HSA-69620",
    "Signaling_by_WNT_R-HSA-195721",
    "PIP3_activates_AKT_signaling_R-HSA-1257604",
    "Mitotic_Prophase_R-HSA-68875",
    "HDR_through_Homologous_Recombination_HRR_or_Single_Strand_Annealing_SSA__R-HSA-5693567",
]


def audit(dirname):
    d = CATALOG / dirname
    if not d.is_dir():
        print(f"### {dirname}: MISSING"); return
    stid = {}
    with open(d / "stid_to_uuid_mapping.csv") as f:
        for r in csv.DictReader(f):
            stid[r["uuid"]] = r["stable_id"]
    deg = defaultdict(int)
    is_target = set()
    reg_cat_edges = []  # (source_uuid, edge_type)
    with open(d / "logic_network.csv") as f:
        for r in csv.DictReader(f):
            s, t, et = r["source_id"], r["target_id"], r["edge_type"]
            deg[s] += 1; deg[t] += 1
            is_target.add(t)
            if et in ("regulator", "catalyst"):
                reg_cat_edges.append((s, et))

    uuids_by_stid = defaultdict(list)
    for u, s in stid.items():
        uuids_by_stid[s].append(u)

    # ---- #8: arbitrary reg/cat attachment ----
    multi = {s: us for s, us in uuids_by_stid.items() if len(us) >= 2}
    suspect = 0
    total_multi_regcat = 0
    for src, et in reg_cat_edges:
        s = stid.get(src)
        if s in multi:
            total_multi_regcat += 1
            best = max(multi[s], key=lambda u: deg[u])
            # suspect if the chosen source is not the max-degree occurrence
            # AND a strictly better occurrence exists
            if deg[src] < deg[best]:
                suspect += 1

    # ---- #7: entity both produced and root ----
    both_prod_and_root = 0
    for s, us in uuids_by_stid.items():
        produced = any(u in is_target for u in us)
        has_root = any(u not in is_target for u in us)
        if produced and has_root and len(us) >= 2:
            both_prod_and_root += 1

    print(f"### {dirname.split('_R-HSA')[0]}")
    print(f"  multi-UUID stIds: {len(multi)}")
    print(f"  #8 reg/cat edges from multi-UUID stIds: {total_multi_regcat}; "
          f"attached to NON-max-degree occurrence: {suspect} "
          f"({100*suspect/max(total_multi_regcat,1):.0f}%)")
    print(f"  #7 stIds both produced-somewhere AND root-elsewhere: {both_prod_and_root}")


def main():
    args = sys.argv[1:]
    dirs = EXP_DIRS
    if args:
        dirs = [d for d in EXP_DIRS if any(a.lower() in d.lower() for a in args)]
    for d in dirs:
        audit(d)


if __name__ == "__main__":
    main()
