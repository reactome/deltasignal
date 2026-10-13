"""specs/021 pilot: how much of a held-out phospho axis decryptM can supply.

A case is (drug-target gene knocked down, phosphosite readout) within one
pathway network. It is scorable when decryptM measured that site under a drug
whose primary target is the gene. The measured direction is the curve's
regulation call: down, up, or not (measured and unchanged).

Inputs:
  BUILD        catalog build dir (stid_to_uuid_mapping.csv, containment.csv per pathway)
  PHOSPHO_CSV  Neo4j export: EWAS stId, UniProt, gene, coordinate, PSI-MOD name
  TARGET_CSV   Neo4j export: EWAS stId, gene, for the drug-target genes
  SITES_TSV    decryptm_sites.py output (any number)

Prints coverage per split and writes one row per scorable (gene, site) pair
to OUT_TSV. With CASES_TSV set, it also writes the scored unit of the
pre-registration (research.md, 2026-10-12): one row per (drug, pathway, site),
a multi-target drug being one perturbation of all its targets in the network.

Usage: [CASES_TSV=cases.tsv] decryptm_coverage.py BUILD PHOSPHO_CSV TARGET_CSV OUT_TSV SITES_TSV [SITES_TSV ...]
"""
import collections
import csv
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from holdout_report import TUNING_PATHWAYS  # noqa: E402

RESIDUE = {"O-phospho-L-serine": "S", "O-phospho-L-threonine": "T", "O4'-phospho-L-tyrosine": "Y"}

# Primary targets only. Staurosporine (pan-kinase) and non-kinase drugs are excluded.
HAND_DRUG_TARGETS = {
    "MK2206": ["AKT1", "AKT2", "AKT3"],
    "Pictilisib": ["PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG"],
    "Dactolisib": ["PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "MTOR"],
    "AZD8055": ["MTOR"],
    "Selumetinib": ["MAP2K1", "MAP2K2"], "PD325901": ["MAP2K1", "MAP2K2"],
    "Refametinib": ["MAP2K1", "MAP2K2"],
    "Gefitinib": ["EGFR"], "Afatinib": ["EGFR", "ERBB2", "ERBB4"], "Lapatinib": ["EGFR", "ERBB2"],
    "Trastuzumab": ["ERBB2"], "Pertuzumab": ["ERBB2"],
    "AZD4547": ["FGFR1", "FGFR2", "FGFR3"],
    "Imatinib": ["ABL1", "KIT", "PDGFRA", "PDGFRB"],
    "Dasatinib": ["ABL1", "SRC", "LCK", "LYN", "FYN", "YES1", "KIT"],
    "SHP099": ["PTPN11"],
    "Nintedanib": ["KDR", "FLT1", "FLT4", "FGFR1", "FGFR2", "FGFR3", "PDGFRA", "PDGFRB"],
    "Tideglusib": ["GSK3B"],
}


def load_drug_targets(path):
    """{drug: [genes]} from a bench/drug_targets.py table (specs/051 R4)."""
    out = collections.defaultdict(list)
    with open(path, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            out[r["drug"]].append(r["gene"])
    return dict(out)


# specs/051 R4: measured affinity (Kd_app <= 1 uM, kinobeads) replaces the hand
# table, which a drug the assay cannot see still falls back to. The map is built
# by bench/drug_targets.py from CC BY-SA data and kept beside that data, outside
# this repository (Adam, 2026-10-13). DRUG_TARGETS=hand reproduces specs/021; any
# other value is a path to a drug_targets.py table.
DEFAULT_DRUG_TARGETS = Path.home() / "deltasignal-catalogs/external/klaeger2017_chembl37/drug_targets_1uM.tsv"


def resolve_drug_targets():
    dt = os.environ.get("DRUG_TARGETS", str(DEFAULT_DRUG_TARGETS))
    if dt == "hand":
        return HAND_DRUG_TARGETS
    if not Path(dt).exists():
        raise SystemExit(f"no drug->target map at {dt}: build it with bench/drug_targets.py "
                         "(specs/051 R4), or set DRUG_TARGETS=hand for the specs/021 table")
    return load_drug_targets(dt)


def neo4j_rows(path):
    for line in open(path):
        if line.startswith('"'):
            yield [x.strip().strip('"') for x in line.split(",")]


def drug_of(experiment):
    """ddPTM_<cell line>_<drug>_<time>_<rep>. A name with "+" is rejected here, but
    decryptM's combinations are mostly written without it (`GeftinibAZD4547-1to80`,
    `SelumetinibMK2206-1to2`); those drop out only because they are not keys of
    the drug->target map (code review 2026-10-13)."""
    parts = experiment.split("_")
    if len(parts) < 4:
        return None, None
    drug = parts[2]
    return (drug if "+" not in drug else None), parts[1]


def call(regulations):
    """Several experiments' calls -> one: regulated calls that all agree give that
    direction, disagreeing ones `conflict`, none `not`."""
    called = set(regulations) - {"not"}
    return "conflict" if len(called) > 1 else (called.pop() if called else "not")


def network_sites(pdir, phospho, targets):
    nodes = {r["stable_id"].split("::")[0] for r in csv.DictReader(open(pdir / "stid_to_uuid_mapping.csv"))}
    cont = collections.defaultdict(set)
    for r in csv.DictReader(open(pdir / "containment.csv")):
        cont[r["stable_id"]].add(r["contains_stable_id"])
    inside = set()
    for n in nodes:
        inside |= cont.get(n, set()) | {n}
    sites = {s for e in inside for s in phospho.get(e, ())}
    genes = {targets[e] for e in inside if e in targets}
    return sites, genes


def main():
    DRUG_TARGETS = resolve_drug_targets()
    build, phospho_csv, target_csv, out, *site_files = sys.argv[1:]
    build = Path(build)
    phospho = collections.defaultdict(set)
    for st, acc, gene, coord, mod in neo4j_rows(phospho_csv):
        if mod in RESIDUE:
            phospho[st].add((acc, int(coord), gene))
    targets = {st: g for st, g in neo4j_rows(target_csv)}
    # measured[(uniprot, position)][gene] -> list of (regulation, drug, cell line)
    measured = collections.defaultdict(lambda: collections.defaultdict(list))
    # by_drug[(uniprot, position)][drug] -> list of (regulation, cell line)
    by_drug = collections.defaultdict(lambda: collections.defaultdict(list))
    for f in site_files:
        for r in csv.DictReader(open(f), delimiter="\t"):
            drug, cell = drug_of(r["experiment"])
            for g in DRUG_TARGETS.get(drug, ()):
                measured[(r["uniprot"], int(r["position"]))][g].append((r["regulation"], drug, cell))
            if drug in DRUG_TARGETS:
                by_drug[(r["uniprot"], int(r["position"]))][drug].append((r["regulation"], cell))
    cases_out = os.environ.get("CASES_TSV")
    fc = open(cases_out, "w", newline="") if cases_out else None
    wc = csv.writer(fc, delimiter="\t") if fc else None
    if wc:
        wc.writerow(["pathway_id", "pathway", "split", "drug", "targets", "uniprot", "position",
                     "site_gene", "measured", "cell_lines"])
    case_tally = collections.Counter()
    names = {}
    for r in csv.DictReader(open(build / "pathways.tsv"), delimiter="\t"):
        v = list(r.values())
        names[v[0]] = v[1]
    tally = collections.Counter()
    with open(out, "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow(["pathway", "split", "gene", "uniprot", "position", "site_gene", "measured", "drugs"])
        for pid in sorted(names):
            pdir = build / pid
            if not pdir.is_dir():
                continue
            split = "tuning" if names[pid] in TUNING_PATHWAYS else "held-out"
            sites, genes = network_sites(pdir, phospho, targets)
            for drug in sorted(DRUG_TARGETS):
                hit = sorted(set(DRUG_TARGETS[drug]) & genes)
                if not hit:
                    continue
                for acc, pos, site_gene in sorted(sites):
                    calls = by_drug.get((acc, pos), {}).get(drug)
                    if not calls:
                        continue
                    m = call(c[0] for c in calls)
                    case_tally[split, m] += 1
                    if wc:
                        wc.writerow([pid, names[pid], split, drug, ";".join(hit), acc, pos, site_gene, m,
                                     ";".join(sorted({c[1] for c in calls}))])
            for gene in sorted(genes):
                tally[split, "possible"] += len(sites)
                for acc, pos, site_gene in sorted(sites):
                    calls = measured.get((acc, pos), {}).get(gene)
                    if not calls:
                        continue
                    # several drugs/lines: a regulated call wins only if all regulated calls agree
                    m = call(c[0] for c in calls)
                    tally[split, m] += 1
                    w.writerow([names[pid], split, gene, acc, pos, site_gene, m,
                                ";".join(sorted({f"{d}@{c}" for _, d, c in calls}))])
    for split in ("held-out", "tuning"):
        n = {k: tally[split, k] for k in ("possible", "down", "up", "not", "conflict")}
        scorable = n["down"] + n["up"] + n["not"]
        print(f"{split}: (gene, site) pairs possible {n['possible']}, measured {scorable + n['conflict']} "
              f"(down {n['down']}, up {n['up']}, unchanged {n['not']}, conflicting {n['conflict']})")
        c = {k: case_tally[split, k] for k in ("down", "up", "not", "conflict")}
        print(f"{split}: (drug, pathway, site) cases {sum(c.values())} "
              f"(down {c['down']}, up {c['up']}, unchanged {c['not']}, conflicting {c['conflict']})")
    if fc:
        fc.close()


if __name__ == "__main__":
    main()
