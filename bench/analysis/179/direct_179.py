#!/usr/bin/env python3
"""specs/179: the curated direct-substrate table. Values-blind (networks and identifier tables only).

A (protein, readout) pair is direct when the protein is a leaf of a catalyst of a reaction whose output carries a
phosphate at one of the readout's registered sites that no input or catalyst of the reaction carries. A protein that
is one of the readout's own genes is left out.

Usage: direct_179.py BUILD OUT_TSV
"""
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parents[1] / "173"))
sys.path.insert(0, str(HERE.parents[1] / "174"))
import lun_a_173 as A  # noqa: E402
from trace_174 import READOUT_ACC  # noqa: E402


def site_readouts(readouts):
    """{(accession, position): readout} and {readout: its own genes}."""
    site, own = {}, {}
    for r in readouts:
        own[r["readout"]] = {s.split(":")[0] for s in r["sites"].split(",")}
        for s in r["sites"].split(","):
            g, pos = s.split(":")
            if g in READOUT_ACC:
                site[(READOUT_ACC[g], int(pos))] = r["readout"]
    return site, own


def direct_pairs(build, acc_gene, site, own):
    """{(gene, readout): sorted pathways}."""
    pairs = defaultdict(set)
    for pdir in sorted(p for p in Path(build).iterdir() if (p / "logic_network.csv").exists()):
        ref, mods = {}, {}
        for r in csv.DictReader(open(pdir / "leaf_refs.csv")):
            if r["reference_stid"].startswith("uniprot:"):
                ref[r["leaf_stid"]] = r["reference_stid"].split(":")[1].split("-")[0]
                mods[r["leaf_stid"]] = r["modifications"]
        leaves = {r["uuid"]: [x for x in r["member_leaves"].split("|") if x in ref]
                  for r in csv.DictReader(open(pdir / "nodes.csv"))}
        ins, cats, outs = defaultdict(set), defaultdict(set), defaultdict(set)
        for e in csv.DictReader(open(pdir / "logic_network.csv")):
            if e["edge_type"] == "input":
                ins[e["target_id"]].add(e["source_id"])
            elif e["edge_type"] == "catalyst":
                cats[e["target_id"]].add(e["source_id"])
            elif e["edge_type"] == "output":
                outs[e["source_id"]].add(e["target_id"])

        def sites(nodes):
            return {(ref[x], int(p)) for u in nodes for x in leaves.get(u, ())
                    for p in re.findall(r"phospho[^;]* at (\d+)", mods[x])}
        for rx, cs in cats.items():
            hit = {site[s] for s in sites(outs[rx]) - sites(ins[rx] | cs) if s in site}
            for u in cs:
                for x in leaves.get(u, ()):
                    g = acc_gene.get(ref[x])
                    for r in hit:
                        if g and g not in own[r]:
                            pairs[(g, r)].add(pdir.name)
    return {k: sorted(v) for k, v in pairs.items()}


def main():
    build, out = sys.argv[1:3]
    acc_gene = {r["accession"]: r["gene"] for r in A.read_tsv(A.MAPS / "pois.tsv") if r["accession"]}
    site, own = site_readouts([r for r in A.read_tsv(A.MAPS / "readouts.tsv") if r["kind"] == "phospho"])
    pairs = direct_pairs(build, acc_gene, site, own)
    with open(out, "w", newline="") as f:
        f.write(f"# specs/179: curated direct substrates on the Lun panel, from build {Path(build).name}. No Lun value read.\n")
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "readout", "pathways"])
        for (g, r), pw in sorted(pairs.items()):
            w.writerow([g, r, ",".join(pw)])
    print(f"{len(pairs)} direct pairs, {len({g for g, _ in pairs})} proteins")


if __name__ == "__main__":
    main()
