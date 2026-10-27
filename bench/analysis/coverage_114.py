#!/usr/bin/env python3
"""specs/114 part 1: how many of a dataset's readouts a catalog build holds.

Graph only: reads site identities, never a measured value. A readout is
- in at protein level if a node's `member_leaves` holds a leaf whose
  `leaf_refs.csv` reference is the readout's UniProt accession, and
- in at site level if a node's `member_leaves` holds a human EWAS of that protein
  phosphorylated at that coordinate (PHOSPHO_CSV, the Release 97 Neo4j export).
Gene symbols (korkut) resolve through GENE_CSV, every human ReferenceGeneProduct's
gene names and accession in Release 97; a symbol may name several accessions.
A readout with several sites is in if any of them is.

Readouts:
  decryptm  every distinct (UniProt, position) in all_sites.tsv (columns 1-2 only)
  korkut    bench/analysis/062/readouts.tsv, kind phospho or total (total: protein level only)
  p100      the 96 analytes (pr_uniprot_id, pr_p100_phosphosite)

Usage: coverage_114.py PHOSPHO_CSV GENE_CSV DATASET READOUT_FILE BUILD [BUILD ...]
"""
import collections
import csv
import gzip
import re
import sys
from pathlib import Path


RESIDUE = {"O-phospho-L-serine", "O-phospho-L-threonine", "O4'-phospho-L-tyrosine"}


def acc(x):
    """UniProt accession without isoform suffix."""
    return x.strip().split("-")[0]


def neo4j_rows(path):
    for line in open(path):
        if line.startswith('"'):
            yield [x.strip().strip('"') for x in line.split(",")]


def phospho_sites(path):
    """(accession, coordinate) -> EWAS stIds, and (gene, coordinate) -> the same."""
    by_acc, by_gene = collections.defaultdict(set), collections.defaultdict(set)
    for stid, uniprot, gene, coord, psi in neo4j_rows(path):
        if coord in ("", "null") or psi not in RESIDUE:
            continue
        by_acc[(acc(uniprot), int(coord))].add(stid)
        by_gene[(gene, int(coord))].add(stid)
    return by_acc, by_gene


def build_contents(build):
    """(EWAS stIds held by some node, accessions held by some node) over every pathway of BUILD."""
    leaves, accs = set(), set()
    for net in sorted(Path(build).glob("R-HSA-*")):
        ref = {}
        for r in csv.DictReader(open(net / "leaf_refs.csv")):
            if r["reference_stid"].lower().startswith("uniprot:"):
                ref[r["leaf_stid"]] = acc(r["reference_stid"].split(":", 1)[1])
        for r in csv.DictReader(open(net / "nodes.csv")):
            for x in filter(None, r["member_leaves"].split("|")):
                leaves.add(x)
                if x in ref:
                    accs.add(ref[x])
    return leaves, accs


def decryptm_readouts(path, by_acc):
    seen = {}
    with open(path) as f:
        next(f)
        for line in f:
            u, pos = line.split("\t", 2)[:2]
            seen[(acc(u), int(pos))] = None
    return [(f"{u}:{p}", {u}, by_acc.get((u, p), set())) for u, p in seen]


def korkut_readouts(path, by_gene, gene_acc):
    out = []
    rows = csv.DictReader((l for l in open(path) if not l.startswith("#")), delimiter="\t")
    for r in rows:
        if r["kind"] == "phospho":
            gcs = [(g, int(c)) for g, c in (s.split(":") for s in r["sites"].split(","))]
            out.append((r["readout"], set().union(*(gene_acc.get(g, set()) for g, _ in gcs)),
                        set().union(*(by_gene.get(gc, set()) for gc in gcs))))
        elif r["kind"] == "total":
            out.append((r["readout"], set().union(*(gene_acc.get(g, set()) for g in r["sites"].split(","))), None))
    return out


def p100_readouts(path, by_acc):
    out = []
    with gzip.open(path, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            u = acc(r["pr_uniprot_id"])
            coords = [int(c) for c in re.findall(r"[STY](\d+)", r["pr_p100_phosphosite"])]
            out.append((r["pr_analyte_id"], {u}, set().union(*(by_acc.get((u, c), set()) for c in coords))))
    return out


def counts(readouts, leaves, accs):
    """(readouts, protein-level in, with a site in Reactome, site-level in)."""
    n = len(readouts)
    prot = sum(bool(a & accs) for _, a, _ in readouts)
    with_site = [s for _, _, s in readouts if s]
    site = sum(bool(s & leaves) for s in with_site)
    return n, prot, len(with_site), site


def main():
    phospho, genes, dataset, readout_file, builds = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5:]
    by_acc, by_gene = phospho_sites(phospho)
    gene_acc = collections.defaultdict(set)
    for g, u in neo4j_rows(genes):
        gene_acc[g].add(acc(u))
    readouts = {"decryptm": lambda: decryptm_readouts(readout_file, by_acc),
                "korkut": lambda: korkut_readouts(readout_file, by_gene, gene_acc),
                "p100": lambda: p100_readouts(readout_file, by_acc)}[dataset]()
    print("dataset\tbuild\treadouts\tprotein_in\tsite_in_reactome\tsite_in")
    for b in builds:
        leaves, accs = build_contents(b)
        print(dataset, Path(b).name, *counts(readouts, leaves, accs), sep="\t")


if __name__ == "__main__":
    main()
