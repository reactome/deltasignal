#!/usr/bin/env python3
"""specs/157: one ensemble variant's calls over the recorded benchmark.

Replays the recorded cases through the SERVED stack on a replay server whose
env carries the canonical serving config PLUS exactly one registry variant
(or none: --variant canonical, the reference run), and writes every readout's
fold and call. The scorer (score_157.py) labels robust/fragile per the
registered R2 rule and runs the gate.

The env is verified POSITIVELY from the solve response's own diagnostics
(the specs/156 container-env rule) — and not just the variant's: every OTHER
registry env must read CANONICAL, so a leftover override can never
contaminate a run.

Axes:
- curator + experimental: the tracked specs/152 g4 dumps (pins, rails,
  readouts), the specs/154 replay pattern — solve per (pathway, gene,
  direction), value the recorded key_outputs via /api/readouts.
- korkut: the recorded cases TSV's (pathway, drug set) solves with the
  recorded pin protocol (drug_activity machinery), readouts via
  /api/readouts with each readout's stable ids (mode per kind).

Usage:
  DELTASIGNAL_BASE=http://localhost:PORT DS_CATALOG_ROOT=<served build> \
    ensemble_157.py --variant NAME REGISTRY_TSV G4_DIR KORKUT_CASES_TSV OUT_DIR
"""
import argparse
import csv
import json
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
import benchmark_vs_mpbiopath as bm  # noqa: E402
import drug_activity as da  # noqa: E402
import dataset_roles  # noqa: E402
import korkut_bench as kb  # noqa: E402


def post(ep, body):
    rq = urllib.request.Request(f"{bm.DS_BASE}/api/{ep}", json.dumps(body).encode(),
                                {"Content-Type": "application/json"})
    return json.loads(urllib.request.urlopen(rq, timeout=600).read())


# env -> (response key, canonical check, variant check given the value).
# Checks are prefix tests on the response's own diagnostic — positive
# verification, never the env.
def _pfx(expected):
    return lambda v: isinstance(v, str) and v.startswith(expected)


def _not_pfx(expected):
    return lambda v: isinstance(v, str) and not v.startswith(expected)


FINGERPRINTS = {
    "DS_DEPLETION_OWN_PRODUCT": ("depletion_own_product", _pfx("suppress_only_source"), _pfx),
    "DS_DEPLETION_COMBINE": ("depletion_combine", _pfx("product"), _pfx),
    "DS_SET_POOL_MODE": ("set_pool_mode", _pfx("product"), _pfx),
    "DS_UBIQUITIN_MODE": ("ubiquitin_rule", _pfx("propagate"), _pfx),
    "DS_POOL_SUPPLY_AND": ("pool_supply_and", _pfx("off"), _pfx),
    "DS_POOL_DRIVE_TWINS": ("pool_drive_twins_rule", _pfx("off"), _pfx),
    "DS_POOL_BASE_WEIGHT": ("pool_base_weight", _pfx("split"), _pfx),
    "DS_CARRIER_LOOP_MOVES": ("carrier_loop_moves", _pfx("drop"), _pfx),
    "DS_CARRIER_LOOP_OPEN": ("carrier_loop_open_rule", _pfx("fed_single_ub"), _pfx),
    "DS_CARRIER_COMPONENT_TWIN": ("carrier_component_rule", _not_pfx("off"), _pfx),
    "DS_METABOLITE_POOLS": ("metabolite_pools_rule", _pfx("off"), _pfx),
    "DS_SELF_FED_MODE": ("self_fed_rule", _pfx("off"), _pfx),
    "DS_SELF_ACTIVATOR_LEAVES": ("self_activator_leaves", _pfx("off"),
                                 lambda v: _not_pfx("off") if v == "1" else _pfx(v)),
    "DS_SELF_ACTIVATOR_WEIGHT": ("self_activator_rule", _not_pfx("off"), _pfx),
    "DS_SELF_INHIBITOR_COPIES": ("self_inhibitor_copies", _pfx("on"), _pfx),
    "DS_CYCLE_REEVAL_CARRIERS": ("cycle_reeval_carrier_rule", _pfx("carrier"), _pfx),
    "DS_RNA_TWIN": ("rna_twin_rule", _pfx("on"), _pfx),
    "DS_DRUG_MODE": ("drug_rule", _pfx("inert"), _pfx),
    "DS_PATHOGEN_MODE": ("pathogen_rule", _pfx("inert"), _pfx),
}
# fed_single_ub startswith fed_single, so the fed_single variant's check must
# be exact-with-boundary rather than a bare prefix:
FED_SINGLE = lambda v: isinstance(v, str) and (v == "fed_single" or v.startswith("fed_single:") or v.startswith("fed_single "))  # noqa: E731


def load_registry(path):
    return {r["variant"]: r for r in csv.DictReader(open(path), delimiter="\t")}


def verify_env(registry, variant_row):
    """A probe solve; every registry env read back from the response."""
    build_root = bm.CATALOG_ROOT.resolve()
    pw = sorted(p.name for p in build_root.iterdir()
                if (p / "stid_to_uuid_mapping.csv").exists())[0]
    net = post("parse", {"pathway_id": pw})
    assert net.get("status") == "success", f"probe parse failed on {pw}"
    u = net["nodes"][0]["uuid"]
    sol = post("solve", {"network_id": net["network_id"], "observations": {u: [50.0, 1.0]}})
    for env, (key, canon_check, var_check) in FINGERPRINTS.items():
        got = sol.get(key)
        if variant_row is not None and variant_row["env"] == env:
            value = variant_row["value"]
            if env == "DS_CARRIER_LOOP_OPEN" and value == "fed_single":
                ok = FED_SINGLE(got)
            elif env == "DS_SELF_ACTIVATOR_LEAVES":
                ok = var_check(value)(got)
            else:
                ok = var_check(value)(got)
            if not ok:
                raise SystemExit(f"env probe FAILED: {env}={value} but {key}={got!r}")
        else:
            if not canon_check(got):
                raise SystemExit(f"env probe FAILED: {key}={got!r} is not canonical "
                                 f"(a leftover {env} override?)")
    print("env probe: variant live, every other rule canonical", flush=True)


def curator_leg(g4, writer):
    pins = {p["pathway_name"]: p for p in json.load(open(g4 / "bench_pins.json"))}
    rails = json.load(open(g4 / "bench_rails.json"))
    reads = json.load(open(g4 / "bench_readouts.json"))
    stids_of = defaultdict(dict)
    for p in reads["curator"]:
        for q in p["readouts"]:
            stids_of[p["pathway_name"]][q["ko"]] = (list(q["stable_ids"]), list(q["members"]))
    cases = defaultdict(set)
    for c in rails["cases"]:
        cases[(c["pathway"], c["gene"], c["direction"])].add((c["axis"], c["key_output"]))
    by_pw = defaultdict(list)
    for (pw, gene, d), kos in cases.items():
        by_pw[pw].append((gene, d, kos))
    nid_cache, solves = {}, 0
    for pw in sorted(by_pw):
        pq = pins.get(pw)
        if pq is None:
            continue
        if pw not in nid_cache:
            nid_cache[pw] = post("parse", {"pathway_id": pq["dir"]})["network_id"]
        nid = nid_cache[pw]
        pin_of = {q["gene"]: q["pins"] for q in pq["queries"]}
        for gene, d, kos in sorted(by_pw[pw]):
            pu = pin_of.get(gene, [])
            if not pu:
                continue
            act = 0.0 if d == "0" else 80.0
            sol = post("solve", {"network_id": nid,
                                 "observations": {u: [act, 0.9] for u in pu}})
            solves += 1
            queries = []
            for axis, ko in sorted(kos):
                st, mem = stids_of[pw].get(ko, ([f"R-HSA-{ko}"], []))
                queries.append({"id": f"{axis}|{ko}", "stable_ids": st, "members": mem})
            ro = post("readouts", {"network_id": nid,
                                   "activities": sol.get("node_activities") or {},
                                   "pinned": pu, "queries": queries,
                                   "converged": bool(sol.get("converged", False))})
            down = ro["protocol"]["down_cutoff"]; up = ro["protocol"]["up_cutoff"]
            for e in ro["results"]:
                axis, ko = e["id"].split("|", 1)
                fold = e.get("fold")
                call = ("" if fold is None else
                        "down" if fold < down else "up" if fold >= up else "not")
                writer.writerow([axis, pw, gene, d, ko, fold, call,
                                 int(bool(sol.get("converged")))])
        print(f"  {pw}: done ({solves} solves)", flush=True)
    return solves


def korkut_leg(cases_tsv, writer):
    dataset_roles.require_role("korkut2015", "fit")
    recorded = set()
    for r in csv.DictReader(open(cases_tsv), delimiter="\t"):
        if r["status"] == "scored":
            recorded.add((r["pathway_id"], r["drug_set"], r["readout"]))
    drugs = {d["code"]: d["genes"].split(",") for d in kb.read_tsv(kb.MAPS / "drugs.tsv")}
    readouts = [r for r in kb.read_tsv(kb.MAPS / "readouts.tsv") if r["kind"] != "exclude"]
    tg = sorted({g for gs in drugs.values() for g in gs})
    rg = sorted({s.split(":")[0] for r in readouts for s in r["sites"].split(",")})
    g2s = bm.neo4j_gene_to_stids(sorted(set(tg) | set(rg)))
    sites = sorted({(s.split(":")[0], int(s.split(":")[1]))
                    for r in readouts if r["kind"] == "phospho" for s in r["sites"].split(",")})
    s2s = kb.neo4j_site_stids(sites)
    r_stids = kb.readout_stids(readouts, kb.total_gene_stids(readouts, g2s), s2s)
    kinds = {r["readout"]: r["kind"] for r in readouts}
    units = da.active_units() if da.needs_units() else {}
    actions = da.drug_actions("korkut", sorted(drugs))
    design = kb.conditions()
    drug_sets = sorted({ds for _, ds in design})
    build_root = bm.CATALOG_ROOT.resolve()
    pathways = sorted(p for p in build_root.iterdir()
                      if (p / "stid_to_uuid_mapping.csv").exists())
    solves = 0
    for pdir in pathways:
        s2u = bm.load_stid_to_uuids(pdir)
        wanted = {name for name in r_stids
                  if any((pdir.name, "+".join(ds), name) in recorded for ds in drug_sets)}
        if not wanted:
            continue
        parsed = None
        for ds in drug_sets:
            ds_name = "+".join(ds)
            names = sorted(name for name in r_stids
                           if (pdir.name, ds_name, name) in recorded)
            if not names:
                continue
            genes = sorted({g for c in ds for g in drugs[c]})
            if parsed is None:
                parsed = bm.parse_pathway_via_ds_api(pdir.name)
                if parsed.get("status") != "success":
                    raise SystemExit(f"{pdir.name}: parse failed")
            action = da.gene_actions([(actions[c], drugs[c]) for c in ds])
            pins, obs, cut, roots = da.perturb(bm, pdir, parsed, genes, g2s,
                                               s2u, units, action)
            if not roots:
                continue
            payload = {"network_id": parsed["network_id"], "observations": obs}
            if cut:
                payload["edge_redirects"] = [list(e) for e in sorted(cut)]
            sol = post("solve", payload)
            if sol.get("status") != "success":
                raise SystemExit(f"{pdir.name} {ds_name}: solve failed")
            solves += 1
            queries = [{"id": name, "stable_ids": sorted(r_stids[name]),
                        "mode": "abundance" if kinds.get(name) == "total" else "forms"}
                       for name in names]
            ro = post("readouts", {"network_id": parsed["network_id"],
                                   "activities": sol.get("node_activities") or {},
                                   "pinned": sorted(obs), "queries": queries,
                                   "converged": bool(sol.get("converged", False))})
            down = ro["protocol"]["down_cutoff"]; up = ro["protocol"]["up_cutoff"]
            for e in ro["results"]:
                fold = e.get("fold")
                call = ("" if fold is None else
                        "down" if fold < down else "up" if fold >= up else "not")
                writer.writerow(["korkut", pdir.name, ds_name, "", e["id"], fold, call,
                                 int(bool(sol.get("converged")))])
        print(f"  {pdir.name}: done ({solves} korkut solves)", flush=True)
    return solves


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", required=True)
    ap.add_argument("registry"); ap.add_argument("g4"); ap.add_argument("korkut_cases")
    ap.add_argument("out_dir")
    ap.add_argument("--axes", default="both", choices=["curator", "korkut", "both"])
    a = ap.parse_args()
    registry = load_registry(a.registry)
    row = None if a.variant == "canonical" else registry[a.variant]
    out_dir = Path(a.out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    served = (json.loads(urllib.request.urlopen(f"{bm.DS_BASE}/api/health", timeout=30).read())
              .get("catalog") or {})
    print(f"server {bm.DS_BASE} serving {served.get('build_id')}", flush=True)
    verify_env(registry, row)

    out = out_dir / f"calls_{a.variant}.tsv"
    with open(out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["axis", "pathway", "case_a", "case_b", "readout", "fold", "call",
                    "converged"])
        n = 0
        if a.axes in ("curator", "both"):
            n += curator_leg(Path(a.g4), w)
        if a.axes in ("korkut", "both"):
            n += korkut_leg(a.korkut_cases, w)
    print(f"{a.variant}: {n} solves -> {out}")


if __name__ == "__main__":
    main()
