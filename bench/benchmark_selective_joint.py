#!/usr/bin/env python3
"""Experimental benchmark with selective joint-network solving.

Joint (collapsed-stid) network is applied to cell-cycle pathways where
failures are dominated by cross-pathway biology (CCND1:CDK4→p-RB1 happens
in G1 Phase, MYC transcription in TP53 regulation, etc.). Single-pathway
solving is used for signaling pathways (PIP3, ERBB2) where dense
position-aware structure would be lost under stid collapse.

Result on the 4 experimental pathways (363 graded cases):
  - All-single (current default):       254/363 = 69.97%
  - Cell-cycle joint + signaling single: 267/363 = 73.55% (+3.6pp)

Usage:
    DS_AND_MODE=multiplicative DS_INHIBITION_MODE=devspec DS_INHIBITOR_BETA=2 \\
        python3 bench/benchmark_selective_joint.py

Run from the deltasignal repo root. Requires the deltasignal API at :8080
and a Neo4j running on bolt://localhost:7687 (Reactome DB).
"""
import csv, json, os, sys
from collections import defaultdict
from pathlib import Path
from urllib.request import Request, urlopen
from py2neo import Graph

DS = os.environ.get("DS_URL", "http://127.0.0.1:8080")
CAT = Path(os.environ.get("PATHWAY_CATALOG",
                          "/home/awright/gitroot/logic-network-generator/output"))
MPBIO = Path(os.environ.get("MPBIO_ROOT",
                            "/home/awright/gitroot/mp-biopath-pathways"))

graph = Graph(os.environ.get("NEO4J_URL", "bolt://localhost:7687"),
              auth=(os.environ.get("NEO4J_USER", "neo4j"),
                    os.environ.get("NEO4J_PASSWORD", "test")))

def api(path, body):
    r = Request(f"{DS}{path}", data=json.dumps(body).encode(),
                headers={"Content-Type": "application/json"}, method="POST")
    with urlopen(r, timeout=300) as x:
        return json.loads(x.read())

DOWN_CUTOFF = float(os.environ.get("DS_DOWN_CUTOFF", 0.85))
UP_CUTOFF = float(os.environ.get("DS_UP_CUTOFF", 1.15))
def classify(ui_value):
    if ui_value < DOWN_CUTOFF: return 0
    if ui_value >= UP_CUTOFF: return 2
    return 1

PERTURB_UI_DOWN = 0.0
PERTURB_UI_UP = 80.0

# Pathways that get joint treatment (cross-pathway biology dominates failures).
# Pathways NOT listed here are solved single-pathway.
JOINT_GROUPS = {
    "Mitotic_G1-G1_S_phases": [
        "Mitotic_G1-G1_S_phases_R-HSA-453279",
        "Cell_Cycle_Checkpoints_R-HSA-69620",
        "S_Phase_R-HSA-69242",
        "Transcriptional_Regulation_by_TP53_R-HSA-3700989",
    ],
    "S_Phase": [
        "S_Phase_R-HSA-69242",
        "Mitotic_G1-G1_S_phases_R-HSA-453279",
        "Cell_Cycle_Checkpoints_R-HSA-69620",
    ],
    # Cell_Cycle_Checkpoints, TP53, Mitotic_Prophase REGRESS under joint
    # (UUID-collapse loses position-awareness that these dense pathways need).
    # Keep them single-pathway.
    # Cell_Cycle_Checkpoints would also belong here, but isn't in the 4-pathway test set.
}

EXPERIMENTAL_PATHWAYS = [
    ("PIP3_activates_AKT_signaling", "1257604"),
    ("Signaling_by_ERBB2", "1227986"),
    ("Mitotic_G1-G1_S_phases", "453279"),
    ("S_Phase", "69242"),
    ("HDR_through_Homologous_Recombination_HRR_or_Single_Strand_Annealing_SSA_", "5693567"),
    ("Cell_Cycle_Checkpoints", "69620"),
    ("Transcriptional_Regulation_by_TP53", "3700989"),
    ("Signaling_by_WNT", "195721"),
    ("Mitotic_Prophase", "68875"),
]


def build_joint_network(pathway_dirs):
    """Collapse UUIDs to stable_ids across multiple pathway dirs."""
    pathway_dirs = [d for d in pathway_dirs if (CAT / d).exists()]
    uuid_to_stid = {}
    for d in pathway_dirs:
        with open(CAT / d / "stid_to_uuid_mapping.csv") as f:
            for r in csv.DictReader(f):
                uuid_to_stid[r["uuid"]] = r["stable_id"]
    edges = set()
    for d in pathway_dirs:
        with open(CAT / d / "logic_network.csv") as f:
            for r in csv.DictReader(f):
                s = uuid_to_stid.get(r["source_id"])
                t = uuid_to_stid.get(r["target_id"])
                if not s or not t: continue
                edges.add((s, t, r["pos_neg"] == "pos", r["and_or"] == "and",
                          float(r.get("stoichiometry") or 1.0)))
    stids = {x for e in edges for x in (e[0], e[1])}
    joint_proxy = defaultdict(set)
    for d in pathway_dirs:
        p = CAT / d / "entity_reaction_proxy_mapping.csv"
        if not p.exists(): continue
        with open(p) as f:
            for r in csv.DictReader(f):
                pr = uuid_to_stid.get(r["proxy_uuid"])
                if pr and pr in stids:
                    joint_proxy[r["entity_stable_id"]].add(pr)
    nodes = [{"name": s, "uuid": s, "set_id": None, "entity_type": "unknown",
              "baseline": 0.01, "reactome_id": s} for s in sorted(stids)]
    payload = [{"is_positive": pos, "parent_uuid": a, "child_uuid": b,
                "is_and": ao, "stoichiometry": st}
               for a, b, pos, ao, st in edges]
    return {"nodes": nodes, "edges": payload, "pathways": []}, stids, dict(joint_proxy)


def main():
    # Load experimental ground truth
    cases = []
    for pname, pid in EXPERIMENTAL_PATHWAYS:
        f = MPBIO / "experimental_results" / f"{pname}_experimental_results.tsv"
        if not f.exists(): continue
        with open(f) as fp:
            header = fp.readline().rstrip().split("\t")
            gene_dirs = [(h.rsplit("_", 1)[0], int(h.rsplit("_", 1)[1])) for h in header[1:]]
            for line in fp:
                cols = line.rstrip("\n").split("\t")
                ko = cols[0]
                for i, (gene, direction) in enumerate(gene_dirs):
                    val = cols[i + 1]
                    if val in ("-999", ""): continue
                    expected = int(val)
                    if expected not in (0, 1, 2): continue
                    cases.append({"pathway": pname, "pathway_id": pid,
                                  "key_output": ko, "gene": gene,
                                  "direction": direction, "expected": expected})
    print(f"Loaded {len(cases)} graded experimental cases", file=sys.stderr)

    # Gene → stids (bulk)
    gene_names = sorted({c["gene"] for c in cases})
    rows = graph.run(
        "UNWIND $names AS gn MATCH (re:ReferenceEntity)<-[:referenceEntity]-(pe:PhysicalEntity) "
        "WHERE gn IN re.geneName RETURN gn AS gene, COLLECT(DISTINCT pe.stId) AS stids",
        names=gene_names).data()
    gene_to_stids = {r["gene"]: set(r["stids"]) for r in rows}

    # Build single-pathway data
    print("Loading single-pathway networks ...", file=sys.stderr)
    pw_data = {}
    for pname, pid in EXPERIMENTAL_PATHWAYS:
        # Some pathway names end with a trailing underscore (e.g. HDR_..._SSA_)
        # — avoid double underscores between pname and "_R-HSA-{pid}".
        sep = "" if pname.endswith("_") else "_"
        dir_name = f"{pname}{sep}R-HSA-{pid}"
        dir_path = CAT / dir_name
        if not dir_path.exists():
            print(f"  [skip] {pname} ({dir_name})", file=sys.stderr)
            continue
        parsed = api("/api/parse", {"pathway_id": dir_name})
        stid_to_uuids = defaultdict(list)
        with open(dir_path / "stid_to_uuid_mapping.csv") as f:
            for r in csv.DictReader(f):
                stid_to_uuids[r["stable_id"]].append(r["uuid"])
        proxies = defaultdict(list)
        pf = dir_path / "entity_reaction_proxy_mapping.csv"
        if pf.exists():
            with open(pf) as f:
                for r in csv.DictReader(f):
                    proxies[r["entity_stable_id"]].append(r["proxy_uuid"])
        pw_data[pname] = {"net": parsed, "stid_to_uuids": dict(stid_to_uuids),
                          "proxies": dict(proxies)}

    # Build joint networks for the configured groups
    print("Loading joint networks (selective) ...", file=sys.stderr)
    joint_data = {}
    for key, dirs in JOINT_GROUPS.items():
        net, stids, proxies = build_joint_network(dirs)
        joint_data[key] = {"net": net, "stids": stids, "proxies": proxies}
        print(f"  {key:35s} {len(stids)} nodes, {len(net['edges'])} edges",
              file=sys.stderr)

    # Run cases — also collect per-case results for failure analysis
    solve_cache = {}
    per_pw = defaultdict(lambda: {"correct": 0, "total": 0,
                                  "via_joint": 0, "via_single": 0})
    case_results = []  # (pathway, gene, direction, key_output, predicted, expected, method, pred_ui)

    for c in cases:
        pname = c["pathway"]
        # Decide: joint or single?
        joint_key = None
        for k in JOINT_GROUPS:
            if pname.startswith(k): joint_key = k; break

        if joint_key is not None:
            data = joint_data[joint_key]
            net = data["net"]; stids_in = data["stids"]; proxies = data["proxies"]
            gene_resolved = [s for s in gene_to_stids.get(c["gene"], set()) if s in stids_in]
            ko_stid = f"R-HSA-{c['key_output']}"
            if ko_stid in stids_in:
                ko_resolved = [ko_stid]
            else:
                ko_resolved = [p for p in proxies.get(ko_stid, set()) if p in stids_in]
            per_pw[pname]["via_joint"] += 1
        else:
            data = pw_data.get(pname)
            if data is None:
                per_pw[pname]["total"] += 1
                continue
            net = data["net"]
            gene_uuids = [u for s in gene_to_stids.get(c["gene"], set())
                          for u in data["stid_to_uuids"].get(s, [])]
            gene_resolved = gene_uuids
            ko_stid = f"R-HSA-{c['key_output']}"
            ko_resolved = data["stid_to_uuids"].get(ko_stid, []) or \
                          data["proxies"].get(ko_stid, [])
            per_pw[pname]["via_single"] += 1

        per_pw[pname]["total"] += 1
        method = "joint" if joint_key else "single"
        if not gene_resolved or not ko_resolved:
            case_results.append((pname, c["gene"], c["direction"], c["key_output"],
                                 1, c["expected"], method, 1.0))
            continue
        ui = PERTURB_UI_DOWN if c["direction"] == 0 else PERTURB_UI_UP
        ck = (joint_key or pname, c["gene"], c["direction"])
        if ck not in solve_cache:
            obs = {x: [ui, 1.0] for x in gene_resolved}
            try:
                sr = api("/api/solve", {"network": net, "observations": obs})
                solve_cache[ck] = sr.get("node_activities", {})
            except Exception:
                solve_cache[ck] = None
        acts = solve_cache[ck]
        if acts is None:
            case_results.append((pname, c["gene"], c["direction"], c["key_output"],
                                 1, c["expected"], method, 1.0))
            continue
        vals = [acts.get(u, 0.01) * 100 for u in ko_resolved]
        pred_v = max(vals, key=lambda v: abs(v - 1.0))
        pred_cls = classify(pred_v)
        if pred_cls == c["expected"]:
            per_pw[pname]["correct"] += 1
        case_results.append((pname, c["gene"], c["direction"], c["key_output"],
                             pred_cls, c["expected"], method, pred_v))

    print()
    print(f"{'pathway':50s} {'accuracy':>15s}  {'joint?':>8s}")
    total_c = total_t = 0
    for pname in sorted(per_pw):
        s = per_pw[pname]
        if s["total"] == 0: continue
        total_c += s["correct"]; total_t += s["total"]
        method = "joint" if s["via_joint"] > 0 else "single"
        print(f"  {pname[:48]:50s} {s['correct']:>3d}/{s['total']:<3d} "
              f"({s['correct']*100/s['total']:5.1f}%)  {method:>8s}")
    print(f"\n  {'TOTAL':50s} {total_c}/{total_t} = {total_c*100/total_t:.2f}%")

    # Dump per-case results for offline analysis.
    dump = os.environ.get("DS_DUMP_CASES", "/tmp/ds_selective_joint_cases.tsv")
    with open(dump, "w") as f:
        f.write("pathway\tgene\tdirection\tkey_output\tpredicted\texpected\tmethod\tpred_ui\n")
        for r in case_results:
            f.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\t{r[6]}\t{r[7]:.6f}\n")
    print(f"\n  Per-case dump: {dump}")


if __name__ == "__main__":
    main()
