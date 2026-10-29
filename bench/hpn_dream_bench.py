#!/usr/bin/env python3
"""specs/112: the frozen model on HPN-DREAM sub-challenge 1A (Hill et al. 2016).

The protocol is pre-registered in specs/112-hpn-dream/spec.md. In short:

- The task: for each (cell line, stimulus) context, score every phosphoprotein
  readout for being a descendant of mTOR, i.e. changed by the mTOR inhibitor
  AZD8055. The truth vectors are the challenge's (`SC1A_TrueDescVectors`).
- `predict` (no truth is read): in every network of the catalog (the pathways
  that act on mTOR, specs/112-hpn-dream/mtor_pathways.tsv), condition A pins the
  stimulus ligands' root forms at DS_PERTURB_UI_UP (none for PBS; all six for
  Serum); condition B adds every root form of MTOR at DS_PERTURB_UI_DOWN. A
  readout is the mean over the unpinned nodes that are, or contain, an EWAS
  phosphorylated at one of its sites (specs/112-hpn-dream/site_map.tsv). Its
  score is |ln(B / A)|, the largest over the networks that hold it, and 0 where
  no network does. The reachability baseline scores 1 where a readout node is
  reachable from an MTOR pin. The model has no cell-line context, so the four
  cell lines get the same scores.
- `score`: AUROC per context against the truth vector (ties count half, as the
  challenge's ROC, which adds a point only where the score changes), the mean
  over the challenge's 30 official contexts, for the model and the baselines.

Usage:
  DS_CATALOG_ROOT=<hpn build> hpn_dream_bench.py predict HPN_DIR PHOSPHO_CSV OUT_DIR
  hpn_dream_bench.py score HPN_DIR OUT_DIR
  hpn_dream_bench.py rank TEAM_AUCS_JSON LEADERBOARD_CSV OUT_DIR

`rank` (amendment 2) reads no truth: the teams' per-context AUROCs from DREAMTools
(dream8/D8C1/data/SC1A_aucs.json), the official SC1A sheet as CSV, and OUT_DIR/score.json.

HPN_DIR is the unpacked challenge download (x_experimental/, x_SC1A_TrueDescVectors/).
PHOSPHO_CSV is the Neo4j export of phosphorylated EWAS (stId, UniProt, gene, coordinate, PSI-MOD name).
"""
import collections
import csv
import json
import math
import os
import sys
from pathlib import Path
from urllib.request import urlopen

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import benchmark_vs_mpbiopath as bm  # noqa: E402
import drug_activity as da  # noqa: E402
import dataset_roles  # noqa: E402
from decryptm_bench import git_state, readout_uuids  # noqa: E402
from decryptm_coverage import RESIDUE, neo4j_rows  # noqa: E402

SPEC = Path(__file__).resolve().parents[1] / "specs" / "112-hpn-dream"
CELL_LINES = ("BT20", "BT549", "MCF7", "UACC812")
STIMULI = {"EGF": ["EGF"], "FGF1": ["FGF1"], "HGF": ["HGF"], "IGF1": ["IGF1"], "Insulin": ["INS"],
           "NRG1": ["NRG1"], "PBS": [], "Serum": ["EGF", "FGF1", "HGF", "IGF1", "INS", "NRG1"]}
TARGET = "MTOR"
# The challenge's final scoring (DREAMTools D8C1, HPNScoringNetwork_ranking.valid_ligands_final):
# BT20 Insulin and BT549 NRG1 had no significant team and were left out.
OFFICIAL = {"UACC812": ["PBS", "Serum", "NRG1", "Insulin", "HGF", "EGF", "FGF1", "IGF1"],
            "BT20": ["IGF1", "PBS", "Serum", "NRG1", "HGF", "EGF", "FGF1"],
            "BT549": ["IGF1", "PBS", "Serum", "Insulin", "HGF", "EGF", "FGF1"],
            "MCF7": ["IGF1", "PBS", "Serum", "NRG1", "Insulin", "HGF", "EGF", "FGF1"]}
EPS = 1e-6
# specs/112 amendment 1: an effect below the solver's tolerance is numerical noise
# (1e-15 to 1e-14 where the pin does not reach), not a ranking; it scores 0.
NOISE = 1e-6


def species_order(hpn_dir, cell):
    """The readouts of a cell line in the order of its MIDAS header, which is the
    order of the truth vectors (DREAMTools `load_species`)."""
    path = Path(hpn_dir) / "x_experimental" / "experimental" / "MIDAS" / f"MD_{cell}_main.csv"
    with open(path) as f:
        header = f.readline().strip().split(",")
    return [h[3:] for h in header if h.startswith("DV:")]


def load_site_map(path=SPEC / "site_map.tsv"):
    out = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        out[r["species"]] = [tuple(p.split(":")) for p in r["sites"].split(";") if p]
    return out


def site_ewas(phospho_csv):
    """(gene, residue+coordinate) -> phosphorylated EWAS stIds."""
    out = collections.defaultdict(set)
    for st, _acc, gene, coord, mod in neo4j_rows(phospho_csv):
        if mod in RESIDUE:
            out[(gene, f"{RESIDUE[mod]}{int(coord)}")].add(st)
    return out


def auroc(scores, classes):
    """Area under the ROC curve, ties counted half (Mann-Whitney). Equals the
    challenge's trapezoid over a ROC that adds a point only where the score
    changes. None if a class is empty."""
    pos = [s for s, c in zip(scores, classes) if c == 1]
    neg = [s for s, c in zip(scores, classes) if c == 0]
    if not pos or not neg:
        return None
    wins = sum((p > n) + 0.5 * (p == n) for p in pos for n in neg)
    return wins / (len(pos) * len(neg))


def predict(hpn_dir, phospho_csv, out_dir):
    role = dataset_roles.require_role("hpn_dream", "score")   # specs/051
    repo = Path(__file__).resolve().parents[1]
    sha, dirty = git_state(repo)
    if dirty and not os.environ.get("DS_ALLOW_DIRTY"):
        raise SystemExit(f"src/ or bench/ is uncommitted, so the result could not be traced:\n{dirty}")
    build = bm.CATALOG_ROOT.resolve()
    served = (json.loads(urlopen(f"{bm.DS_BASE}/api/health", timeout=30).read()).get("catalog") or {})
    if served.get("build_id") != build.name:
        raise SystemExit(f"the API serves {served.get('build_id')!r}, not {build.name!r}")
    species = list(dict.fromkeys(s for c in CELL_LINES for s in species_order(hpn_dir, c)))
    sites = load_site_map()
    ewas = site_ewas(phospho_csv)
    stids = {sp: set().union(*[ewas.get(p, set()) for p in sites.get(sp, [])]) for sp in species}
    genes = sorted({TARGET} | {g for gs in STIMULI.values() for g in gs})
    gene_to_stids = bm.neo4j_gene_to_stids(genes)
    pathways = [r["id"] for r in csv.DictReader(open(build / "pathways.tsv"), delimiter="\t")]
    units = da.active_units() if da.DRUG_ACTION == "activity" else {}
    print(f"Protocol: HPN-DREAM SC1A, DS_DRUG_ACTION={da.DRUG_ACTION} DS_PIN_SCOPE={bm.PIN_SCOPE} "
          f"DS_PIN_ROOT_IGNORE={bm.PIN_ROOT_IGNORE or 'none'} DS_PERTURB_UI_UP={bm.PERTURB_UI_UP:g} "
          f"DS_PERTURB_UI_DOWN={bm.PERTURB_UI_DOWN:g}; build {build.name}; code {sha}; "
          f"dataset hpn_dream registered {role}", flush=True)
    rows, solves = [], 0
    rules = collections.Counter()

    def solve(parsed, pins, pid, payload=None, extra=None):
        nonlocal solves
        obs = {u: [v, bm.PIN_CONFIDENCE] for u, v in sorted(pins.items())}
        obs.update(extra or {})
        res = bm.solve_via_ds_api(payload, obs, network_id=None if payload else parsed["network_id"])
        if res.get("status") != "success":
            raise SystemExit(f"{pid}: solve failed: {res.get('message')}")
        bm.refuse_unless_default_model(res, pid)
        solves += 1
        rules.update(f"{k}={v}" for k, v in res.items() if k.endswith("_rule"))
        return res.get("node_activities", {})

    for pid in pathways:
        pdir = build / pid
        if not (pdir / "logic_network.csv").is_file():
            continue
        stid_to_uuids = bm.load_stid_to_uuids(pdir)
        g2u = {g: [u for s in gene_to_stids.get(g, []) for u in stid_to_uuids.get(s, [])] for g in genes}
        pins, _ = bm.pin_occurrences(pdir, g2u, gene_to_stids)
        target_pins = set(pins.get(TARGET, []))
        if not target_pins:
            print(f"{pid}: no root form of {TARGET}; not used", flush=True)
            continue
        parsed = bm.parse_pathway_via_ds_api(pdir.name)
        if parsed.get("status") != "success":
            raise SystemExit(f"{pid}: parse failed: {parsed.get('message')}")
        reach = bm.reachable_from(bm.build_adjacency(pdir), target_pins)
        for stim, ligands in STIMULI.items():
            a_pins = {u: bm.PERTURB_UI_UP for g in ligands for u in pins.get(g, [])}
            b_pins = {**a_pins, **{u: bm.PERTURB_UI_DOWN for u in target_pins}}
            # A is solved even with no ligand pinned, so B is compared with the solver's own baseline, not a
            # constant: their 1e-15 differences are rounding, not effects (amendment 4)
            act_a = solve(parsed, a_pins, pid)
            if da.DRUG_ACTION == "knockout":
                act_b = solve(parsed, b_pins, pid)
            else:      # specs/122: the inhibitor cuts the target's catalyst / regulator edges; nothing is pinned
                _, cut_obs, payload, _ = da.perturb(bm, pdir, parsed, [TARGET], gene_to_stids, stid_to_uuids, units)
                b_pins = dict(a_pins)
                act_b = solve(parsed, a_pins, pid, payload, cut_obs)
            for sp in species:
                ro = [u for u in readout_uuids(stids[sp], stid_to_uuids) if u not in b_pins]
                if not ro:
                    continue
                va = 100.0 * sum(act_a.get(u, 0.01) for u in ro) / len(ro)
                vb = 100.0 * sum(act_b.get(u, 0.01) for u in ro) / len(ro)
                rows.append({"pathway": pid, "stimulus": stim, "species": sp, "readout_nodes": len(ro),
                             "ligands_pinned": sum(1 for g in ligands if pins.get(g)),
                             "a_ui": f"{va:.6g}", "b_ui": f"{vb:.6g}",
                             "effect": f"{abs(math.log((vb + EPS) / (va + EPS))):.6g}",
                             "reach": int(any(u in reach for u in ro))})
        print(f"{pid}: {len(target_pins)} {TARGET} pins", flush=True)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "hpn_readouts.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]) if rows else ["pathway"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    best = {}
    for r in rows:
        k = (r["stimulus"], r["species"])
        e, rc = float(r["effect"]), r["reach"]
        b = best.get(k, (0.0, 0, ""))
        best[k] = (max(b[0], e), max(b[1], rc), b[2] if b[0] >= e else r["pathway"])
    with open(out / "hpn_predictions.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["stimulus", "species", "ds_score", "reach_score", "best_pathway"])
        for stim in STIMULI:
            for sp in species:
                e, rc, p = best.get((stim, sp), (0.0, 0, ""))
                w.writerow([stim, sp, f"{e:.6g}", rc, p])
    meta = {"build": build.name, "code": sha, "dataset": "hpn_dream", "role": role, "solves": solves,
            "rules": dict(rules), "species": len(species),
            "species_in_networks": len({r["species"] for r in rows})}
    (out / "predict.json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta, indent=1))


def score(hpn_dir, out_dir):
    role = dataset_roles.require_role("hpn_dream", "score")
    out = Path(out_dir)
    pred = {(r["stimulus"], r["species"]): r for r in
            csv.DictReader(open(out / "hpn_predictions.tsv"), delimiter="\t")}
    truth_dir = Path(hpn_dir) / "x_SC1A_TrueDescVectors" / "TrueDescVectors"
    results = []
    for cell in CELL_LINES:
        order = species_order(hpn_dir, cell)
        for stim in STIMULI:
            vec = open(truth_dir / f"TrueVec_{cell}_{stim}.csv").read().strip().split(",")
            if len(vec) != len(order):
                raise SystemExit(f"{cell} {stim}: truth has {len(vec)} entries, header {len(order)}")
            keep = [(sp, int(v)) for sp, v in zip(order, vec) if v != "NaN"]
            classes = [c for _, c in keep]
            r = {"cell": cell, "stimulus": stim, "official": stim in OFFICIAL[cell],
                 "n": len(keep), "positives": sum(classes)}
            for name, col in (("deltasignal", "ds_score"), ("reachability", "reach_score")):
                vals = [float(pred[(stim, sp)][col]) for sp, _ in keep]
                r[name] = auroc([v if v >= NOISE else 0.0 for v in vals], classes)
            r["constant"] = auroc([0.0] * len(keep), classes)
            results.append(r)
    off = [r for r in results if r["official"]]
    mean = lambda rs, k: sum(r[k] for r in rs if r[k] is not None) / sum(r[k] is not None for r in rs)  # noqa: E731
    summary = {"dataset": "hpn_dream", "role": role, "contexts": len(results), "official_contexts": len(off),
               "mean_auroc_official": {k: mean(off, k) for k in ("deltasignal", "reachability", "constant")},
               "mean_auroc_all": {k: mean(results, k) for k in ("deltasignal", "reachability", "constant")},
               "per_context": results}
    (out / "score.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps({k: v for k, v in summary.items() if k != "per_context"}, indent=1))


def rank_min(values):
    """Rank descending, ties at the lowest rank (pandas method="min", as DREAMTools)."""
    return [1 + sum(v > x for v in values) for x in values]


def standings(table, contexts):
    """Mean AUROC and mean rank of every entry over `contexts` (a list, repeats allowed)."""
    names = list(table)
    ranks = {n: 0.0 for n in names}
    for ctx in contexts:
        for n, r in zip(names, rank_min([table[n][ctx] for n in names])):
            ranks[n] += r
    k = len(contexts)
    return ({n: sum(table[n][c] for c in contexts) / k for n in names}, {n: ranks[n] / k for n in names})


def place(scores, name, higher_better):
    """1 + number of entries strictly better than `name`."""
    me = scores[name]
    return 1 + sum((v > me) if higher_better else (v < me) for v in scores.values())


SHEET_ALIASES = {"MLBHKKG": "ML_BHK_KG"}


def rank(team_aucs_json, leaderboard_csv, out_dir, boot=10000, seed=112):
    import numpy as np
    out = Path(out_dir)
    ours = json.loads((out / "score.json").read_text())
    contexts = [f"{c}_{s}" for c in CELL_LINES for s in OFFICIAL[c]]
    raw = json.loads(Path(team_aucs_json).read_text())
    teams = {t: {ctx: raw[ctx][t] for ctx in contexts} for t in raw[contexts[0]]}
    sheet = {SHEET_ALIASES.get(r["Team name"], r["Team name"]): r
             for r in csv.DictReader(open(leaderboard_csv)) if r.get("Team name")}
    # Validation (amendment 3): the per-context AUROCs must reproduce the official sheet, which truncates AUROC
    # to four decimals and whose mean ranks include entries later removed (whole thirtieths, never fewer).
    mean_auc, mean_rank = standings(teams, contexts)
    bad = []
    for t, r in sheet.items():
        if t not in teams:
            bad.append({"team": t, "problem": "not in the per-context file"})
            continue
        a, ra = float(r["mean AUROC"]), float(r["mean Rank"])
        auc_ok = a - 1e-9 <= mean_auc[t] < a + 1e-4 + 1e-9
        rank_ok = any(abs(ra - mean_rank[t] - k / len(contexts)) <= 1e-3 for k in range(4))
        if not (auc_ok and rank_ok):
            bad.append({"team": t, "mean_auroc": [mean_auc[t], a], "mean_rank": [mean_rank[t], ra]})
    on_both = sorted(t for t in teams if t in sheet)
    # every pair, not adjacent ones (ties make an adjacent check non-transitive and dependent on input order); teams
    # tied on our mean rank may stand in either order: the removed entries break such ties on the sheet. Each team in
    # an inverted pair counts against the two-team budget.
    inverted = sorted({t for x in on_both for y in on_both for t in (x, y)
                       if mean_rank[x] < mean_rank[y] - 1e-9 and int(sheet[x]["Final rank"]) > int(sheet[y]["Final rank"])})
    bad += [{"team": t, "problem": "order differs from the sheet's final ranks"} for t in inverted]
    bad += [{"team": t, "problem": "not on the sheet"} for t in teams if t not in sheet]
    result = {"teams": len(teams), "sheet_teams": len(sheet), "contexts": len(contexts), "validation_failures": bad}
    if len(bad) > 2:
        result["note"] = "more than two teams fail validation; no rank reported (amendments 2-3)"
    else:
        per = {(r["cell"], r["stimulus"]): r for r in ours["per_context"]}
        rng = np.random.default_rng(seed)
        draws = [[contexts[i] for i in rng.integers(0, len(contexts), len(contexts))] for _ in range(boot)]
        for entry in ("deltasignal", "reachability"):
            table = dict(teams)
            table[entry] = {f"{c}_{s}": per[(c, s)][entry] for c in CELL_LINES for s in OFFICIAL[c]}
            auc, rk = standings(table, contexts)
            b_auc, b_rank = [], []
            for d in draws:
                a, r = standings(table, d)
                b_auc.append(place(a, entry, True))
                b_rank.append(place(r, entry, False))
            result[entry] = {
                "entries": len(table), "mean_auroc": auc[entry], "mean_rank": rk[entry],
                "place_by_mean_rank": place(rk, entry, False), "place_by_mean_auroc": place(auc, entry, True),
                "boot": boot, "seed": seed,
                "place_by_mean_rank_95": [float(np.percentile(b_rank, 2.5)), float(np.percentile(b_rank, 97.5))],
                "place_by_mean_auroc_95": [float(np.percentile(b_auc, 2.5)), float(np.percentile(b_auc, 97.5))]}
    (out / "rank.json").write_text(json.dumps(result, indent=1))
    print(json.dumps(result, indent=1))


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    if mode == "predict" and len(sys.argv) == 5:
        return predict(*sys.argv[2:5])
    if mode == "score" and len(sys.argv) == 4:
        return score(*sys.argv[2:4])
    if mode == "rank" and len(sys.argv) == 5:
        return rank(*sys.argv[2:5])
    raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
