"""specs/113: LINCS P100 compound -> target proteins, by the rule registered there.

Each compound is found in ChEMBL 37 by its standard InChIKey (then the connectivity block, then name or synonym)
and taken to its parent molecule. Targets come from the first source that gives any:
  1. kinobeads (specs/051 R4: the Klaeger deposit, Kd_app <= cutoff, single human proteins);
  2. ChEMBL drug_mechanism (human targets, complexes and families expanded to their components);
  3. ChEMBL binding activities (assay type B, pChEMBL >= 6, single human proteins).
A compound with none is listed and not scored. No P100 value is read.

Usage: p100_targets.py PERT_INFO_GZ KLAEGER_DIR CACHE_DIR [CUTOFF_NM] > map.tsv

ChEMBL responses are cached as JSON under CACHE_DIR, so a rebuild reads no network. The map is derived from
CC BY-SA data and stays beside it, not in this repository (as drug_targets.py).
"""
import csv
import gzip
import hashlib
import json
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

API = "https://www.ebi.ac.uk/chembl/api/data/"
PCHEMBL_MIN = 6.0


def get(cache, endpoint, params):
    """GET API/endpoint.json?params, cached by a hash of the query; follows pagination."""
    q = urllib.parse.urlencode(sorted(params.items()))
    path = Path(cache) / f"{endpoint.replace('/', '_')}_{hashlib.sha256(q.encode()).hexdigest()[:16]}.json"
    if path.exists():
        return json.loads(path.read_text())
    url, out = f"{API}{endpoint}.json?{q}", None
    while url:
        with urllib.request.urlopen(url, timeout=60) as r:
            page = json.load(r)
        time.sleep(0.2)
        if "page_meta" not in page:
            out = page
            break
        key = next(k for k in page if k != "page_meta")
        out = (out or []) + page[key]
        nxt = page["page_meta"].get("next")
        url = urllib.parse.urljoin(API, nxt) if nxt else None
    path.write_text(json.dumps({"query": f"{endpoint}?{q}", "data": out}))
    return {"query": f"{endpoint}?{q}", "data": out}


def parent(mol):
    h = mol.get("molecule_hierarchy") or {}
    return h.get("parent_chembl_id") or mol["molecule_chembl_id"]


def find_molecules(cache, inchikey, name):
    """[(parent chembl id, how)] in the registered order: exact InChIKey, connectivity block, name, synonym.
    Several hits in one step: the lowest parent id, deterministically."""
    tries = [("inchikey", {"molecule_structures__standard_inchi_key": inchikey, "limit": 100})]
    if inchikey and "-" in inchikey:
        tries.append(("connectivity", {"molecule_structures__standard_inchi_key__startswith": inchikey.split("-")[0],
                                       "limit": 100}))
    tries += [("name", {"pref_name__iexact": name, "limit": 100}),
              ("synonym", {"molecule_synonyms__molecule_synonym__iexact": name, "limit": 100})]
    out = []
    for how, params in tries:
        if not next(iter(params.values())):
            continue
        hits = get(cache, "molecule", params)["data"]
        if hits:
            mol = min((parent(m) for m in hits), key=lambda c: int(c[6:]))
            if mol not in (m for m, _ in out):
                out.append((mol, how))
    return out


def target_components(cache, tid):
    """(target_type, [accessions]) for a human target, else (type, [])."""
    t = get(cache, f"target/{tid}", {})["data"]
    if t.get("organism") != "Homo sapiens":
        return t.get("target_type"), []
    return t.get("target_type"), sorted({c["accession"] for c in t.get("target_components", [])
                                         if c.get("component_type") == "PROTEIN" and c.get("accession")})


def kinobeads(klaeger_dir, cutoff_nm):
    """molecule chembl id -> {accession: Kd} under the R4 rule."""
    tgt = {}
    with open(Path(klaeger_dir) / "targets.tsv", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["target_type"] == "SINGLE PROTEIN" and r["organism"] == "Homo sapiens":
                tgt[r["target_chembl_id"]] = r["accession"]
    out, seen = defaultdict(dict), set()
    with open(Path(klaeger_dir) / "activities.tsv", newline="") as f:
        for a in csv.DictReader(f, delimiter="\t"):
            seen.add(a["molecule_chembl_id"])
            acc = tgt.get(a["target_chembl_id"])
            if acc and a["standard_type"] == "Kd" and a["standard_units"] == "nM" and a["standard_relation"] == "=":
                kd = float(a["standard_value"])
                if kd <= cutoff_nm and kd < out[a["molecule_chembl_id"]].get(acc, float("inf")):
                    out[a["molecule_chembl_id"]][acc] = kd
    return out, seen


def targets_of(cache, mol, kb, kb_mols):
    """({accession: Kd or ""}, source) by the registered order of sources."""
    if mol in kb_mols:
        accs = {a: f"{kd:g}" for a, kd in kb[mol].items()}
        if accs:
            return accs, "kinobeads"
    accs = {}
    for m in get(cache, "mechanism", {"molecule_chembl_id": mol, "limit": 1000})["data"]:
        if m.get("target_chembl_id"):
            accs.update({a: "" for a in target_components(cache, m["target_chembl_id"])[1]})
    if accs:
        return accs, "mechanism"
    acts = get(cache, "activity", {"molecule_chembl_id": mol, "assay_type": "B",
                                   "pchembl_value__gte": PCHEMBL_MIN, "limit": 1000})["data"]
    for tid in sorted({a["target_chembl_id"] for a in acts}):
        ttype, comps = target_components(cache, tid)
        if ttype == "SINGLE PROTEIN":
            accs.update({a: "" for a in comps})
    return (accs, "binding") if accs else ({}, "none")


def build(pert_info_gz, klaeger_dir, cache, cutoff_nm=1000.0):
    Path(cache).mkdir(parents=True, exist_ok=True)
    kb, kb_mols = kinobeads(klaeger_dir, cutoff_nm)
    with gzip.open(pert_info_gz, "rt") as f:
        perts = [r for r in csv.DictReader(f, delimiter="\t") if r["pert_type"] == "trt_cp"]
    rows, log = [], []
    for p in sorted(perts, key=lambda r: r["pert_iname"]):
        # amendment 1: a molecule with no target from any source passes to the next identification step
        found = find_molecules(cache, p["inchi_key"], p["pert_iname"])
        mol, how, accs, src = None, "unmatched", {}, "none"
        for mol, how in found:
            accs, src = targets_of(cache, mol, kb, kb_mols)
            if accs:
                break
        log.append((p["pert_iname"], p["pert_id"], mol or "", how, src, len(accs)))
        rows += [(p["pert_iname"], p["pert_id"], mol, a, accs[a], src) for a in sorted(accs)]
    return rows, log


def main(argv):
    if len(argv) not in (4, 5):
        raise SystemExit(__doc__)
    rows, log = build(argv[1], argv[2], argv[3], float(argv[4]) if len(argv) == 5 else 1000.0)
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    w.writerow(["compound", "pert_id", "chembl_parent", "uniprot", "kd_app_nM", "source"])
    w.writerows(rows)
    for r in log:
        print("\t".join(map(str, r)), file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv)
