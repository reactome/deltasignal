#!/usr/bin/env python3
"""specs/021: score the solver against decryptM drug responses.

The protocol is pre-registered in specs/021-empirical-holdout-axis/research.md
("Pre-registration, 2026-10-12"). In short:

- A case is one (drug, pathway, phosphosite) from
  `bench/analysis/decryptm_coverage.py` (CASES_TSV).
- A drug is one solve per pathway. Every root form of every primary target in
  the network is pinned to DS_PERTURB_UI_DOWN at once, by the curator
  bench's own rule (`pin_occurrences`).
- The readout is the mean over the nodes that are, or contain, an EWAS
  phosphorylated at that (UniProt, coordinate), with pinned nodes left out.
  A case whose readout nodes are all pinned is *trivial*: it is counted, not
  scored.
- It is scored against two baselines: always DOWN, and reachability (DOWN when
  any readout node is reachable from a pin, else NORMAL).

Usage: DS_CATALOG_ROOT=<build> decryptm_bench.py CASES_TSV PHOSPHO_CSV OUT_DIR

PHOSPHO_CSV is the Neo4j export the case table was built from. OUT_DIR gets
decryptm_cases.tsv (one row per case, every status) and summary.json.
"""
import collections
import csv
import json
import math
import os
import subprocess
import sys
from pathlib import Path
from urllib.request import urlopen

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import benchmark_vs_mpbiopath as bm  # noqa: E402
from decryptm_coverage import RESIDUE, neo4j_rows  # noqa: E402

# ERK1/2 activation loop: MAPK3 T202/Y204, MAPK1 T185/Y187 (registered control)
ERK_LOOP = {("P27361", 202), ("P27361", 204), ("P28482", 185), ("P28482", 187)}
NAMES = {bm.DOWN: "down", bm.NORMAL: "not", bm.UP: "up"}
CODE = {v: k for k, v in NAMES.items()}


def readout_uuids(site_stids, stid_to_uuids):
    """Every node that is, or contains, one of the site's phospho EWAS, in order."""
    return list(dict.fromkeys(u for s in sorted(site_stids) for u in stid_to_uuids.get(s, ())))


def score_case(readout, pinned, activities, reach):
    """(status, pred_ui, predicted, reachability prediction) for one case."""
    if not pinned:
        return "no_pin", None, None, None
    if not readout:
        return "no_readout", None, None, None
    free = [u for u in readout if u not in pinned]
    if not free:
        return "trivial", None, None, None
    pred_ui = 100.0 * sum(activities.get(u, 0.01) for u in free) / len(free)
    reach_pred = bm.DOWN if any(u in reach for u in free) else bm.NORMAL
    return "scored", pred_ui, bm.classify(pred_ui), reach_pred


def mcnemar_p(b, c):
    """Exact two-sided McNemar p from the discordant counts."""
    n = b + c
    if n == 0:
        return 1.0
    tail = sum(math.comb(n, i) for i in range(min(b, c) + 1)) / 2 ** n
    return min(1.0, 2 * tail)


def compare(rows):
    """DeltaSignal, reachability and always-DOWN on the same scored cases."""
    n = len(rows)
    ds = [r["predicted"] == r["measured"] for r in rows]
    rc = [r["reach_predicted"] == r["measured"] for r in rows]
    b = sum(x and not y for x, y in zip(ds, rc))
    c = sum(y and not x for x, y in zip(ds, rc))
    return {"n": n, "deltasignal": sum(ds), "reachability": sum(rc),
            "always_down": sum(r["measured"] == "down" for r in rows),
            "ds_only": b, "reach_only": c, "mcnemar_p": mcnemar_p(b, c)}


def macro_f1(rows, key):
    conf = collections.Counter((CODE[r[key]], CODE[r["measured"]]) for r in rows)
    return bm.summary_metrics(conf)["macro_f1"]


def summarize(rows, split="held-out"):
    """The registered metrics for one split."""
    rows = [r for r in rows if r["split"] == split]
    out = {"status": dict(collections.Counter(r["status"] for r in rows))}
    scored = [r for r in rows if r["status"] == "scored"]
    reg = [r for r in scored if r["measured"] in ("down", "up")]
    erk = lambda r: (r["uniprot"], int(r["position"])) in ERK_LOOP  # noqa: E731
    das = lambda r: r["drug"] == "Dasatinib"  # noqa: E731
    out["regulated"] = compare(reg)
    out["regulated_no_erk"] = compare([r for r in reg if not erk(r)])
    out["regulated_no_dasatinib"] = compare([r for r in reg if not das(r)])
    out["regulated_neither"] = compare([r for r in reg if not erk(r) and not das(r)])
    out["up"] = compare([r for r in reg if r["measured"] == "up"])
    out["down"] = compare([r for r in reg if r["measured"] == "down"])
    out["unchanged"] = compare([r for r in scored if r["measured"] == "not"])
    for r in scored:
        r["always_down"] = "down"
    out["macro_f1"] = {k: macro_f1(scored, k) for k in ("predicted", "reach_predicted", "always_down")
                       } if scored else {}
    for key, name in ((lambda r: r["drug"], "by_drug"), (lambda r: r["pathway"], "by_pathway"),
                      (lambda r: f'{r["site_gene"]} {r["uniprot"]} {r["position"]}', "by_site")):
        groups = collections.defaultdict(list)
        for r in reg:
            groups[key(r)].append(r)
        out[name] = {k: compare(v) for k, v in sorted(groups.items())}
    return out


def git_state(repo):
    sha = subprocess.run(["git", "-C", str(repo), "rev-parse", "--short", "HEAD"],
                         capture_output=True, text=True, check=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(repo), "status", "--porcelain", "--", "src", "bench"],
                           capture_output=True, text=True, check=True).stdout.strip()
    return sha, dirty


def main():
    cases_tsv, phospho_csv, out_dir = sys.argv[1:4]
    repo = Path(__file__).resolve().parents[1]
    sha, dirty = git_state(repo)
    if dirty and not os.environ.get("DS_ALLOW_DIRTY"):
        raise SystemExit(f"src/ or bench/ is uncommitted, so the result could not be traced:\n{dirty}")
    build = bm.CATALOG_ROOT.resolve()
    served = (json.loads(urlopen(f"{bm.DS_BASE}/api/health", timeout=30).read()).get("catalog") or {})
    if served.get("build_id") != build.name:
        raise SystemExit(f"the API serves {served.get('build_id')!r}, not {build.name!r}")
    site_stids = collections.defaultdict(set)
    for st, acc, _gene, coord, mod in neo4j_rows(phospho_csv):
        if mod in RESIDUE:
            site_stids[(acc, int(coord))].add(st)
    cases = list(csv.DictReader(open(cases_tsv), delimiter="\t"))
    genes = sorted({g for c in cases for g in c["targets"].split(";")})
    gene_to_stids = bm.neo4j_gene_to_stids(genes)
    by_pathway = collections.defaultdict(lambda: collections.defaultdict(list))
    for c in cases:
        by_pathway[c["pathway_id"]][(c["drug"], c["targets"])].append(c)
    print(f"Protocol: decryptM drug bench, DS_PIN_SCOPE={bm.PIN_SCOPE} "
          f"DS_PERTURB_UI_DOWN={bm.PERTURB_UI_DOWN:g} DS_KO_AGG=mean(unpinned) "
          f"cutoffs {bm.DOWN_CUTOFF}/{bm.UP_CUTOFF}; build {build.name}; code {sha}", flush=True)
    rows, solves = [], 0
    for pid in sorted(by_pathway):
        pdir = build / pid
        stid_to_uuids = bm.load_stid_to_uuids(pdir)
        parsed = bm.parse_pathway_via_ds_api(pdir.name)
        if parsed.get("status") != "success":
            raise SystemExit(f"{pid}: parse failed: {parsed.get('message')}")
        server = {str(n["uuid"]) for n in parsed["nodes"]}
        local = {u for us in stid_to_uuids.values() for u in us}
        if local and not server & local:
            raise SystemExit(f"{pid}: the API's network shares no uuid with {pdir}")
        adj = bm.build_adjacency(pdir)
        for (drug, targets), group in sorted(by_pathway[pid].items()):
            g2u = {g: [u for s in gene_to_stids.get(g, []) for u in stid_to_uuids.get(s, [])]
                   for g in targets.split(";")}
            pins, _ = bm.pin_occurrences(pdir, g2u, gene_to_stids)
            pinned = set().union(*pins.values()) if pins else set()
            activities, reach = {}, set()
            if pinned and any(c["measured"] != "conflict" for c in group):
                obs = {u: [bm.PERTURB_UI_DOWN, bm.PIN_CONFIDENCE] for u in sorted(pinned)}
                res = bm.solve_via_ds_api(None, obs, network_id=parsed["network_id"])
                if res.get("status") != "success":
                    raise SystemExit(f"{pid} {drug}: solve failed: {res.get('message')}")
                bm.refuse_unless_default_model(res, pid)
                solves += 1
                activities = res.get("node_activities", {})
                reach = bm.reachable_from(adj, pinned)
            for c in group:
                ro = readout_uuids(site_stids.get((c["uniprot"], int(c["position"])), ()), stid_to_uuids)
                if c["measured"] == "conflict":
                    status, ui, pred, rp = "conflict", None, None, None
                else:
                    status, ui, pred, rp = score_case(ro, pinned, activities, reach)
                rows.append({**c, "status": status, "pinned": len(pinned), "readout_nodes": len(ro),
                             "pred_ui": "" if ui is None else f"{ui:.6g}",
                             "predicted": NAMES.get(pred, ""), "reach_predicted": NAMES.get(rp, "")})
        print(f"  {pid} {group[0]['pathway']}: {len(by_pathway[pid])} drugs", flush=True)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cols = list(rows[0].keys())
    with open(out / "decryptm_cases.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    summary = {"build": build.name, "code": sha, "solves": solves,
               "held-out": summarize(rows, "held-out"), "tuning": summarize(rows, "tuning")}
    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    for split in ("held-out", "tuning"):
        s = summary[split]
        print(f"{split}: status {s['status']}")
        for k in ("regulated", "regulated_no_erk", "regulated_no_dasatinib", "regulated_neither",
                  "down", "up", "unchanged"):
            m = s[k]
            print(f"  {k:24s} n={m['n']:4d}  DS {m['deltasignal']:4d}  reach {m['reachability']:4d}  "
                  f"always-down {m['always_down']:4d}  (DS-only {m['ds_only']}, reach-only "
                  f"{m['reach_only']}, McNemar p {m['mcnemar_p']:.3g})")
        print(f"  macro-F1 {s['macro_f1']}")


if __name__ == "__main__":
    main()
