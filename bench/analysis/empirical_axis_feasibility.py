#!/usr/bin/env python3
"""How much of the curator axis could a DERIVED empirical axis ever score?

The experimental ground truth in this project covers ten pathways, and all ten
are the tuning ten, so there is no held-out empirical test for either tool.
That is the single biggest threat to the headline claim: we reproduce curator
*reasoning* far better than MP-BioPath and are level with it on experimental
*outcomes*, which is what a model tracking the representation rather than the
biology would look like.

MP-BioPath's experimental truth is a hand-curated 0/1/2 matrix per pathway, so
extending it their way means manual curation. This script sizes the
alternative: deriving an empirical axis from public perturbation data.

It answers one question -- how many (pathway, perturbed gene, readout) triples
in pathways WITHOUT experimental truth have a readout that a transcriptional
assay could actually measure -- and reports the funnel honestly, because most
readouts cannot be measured that way:

  * a readout that is a COMPLEX has no transcript of its own;
  * a phospho- or ubiquitin-form is a state change, not an abundance change;
  * a proteolytic fragment is not a distinct transcript.

Only a readout that is a bare gene symbol survives. That is a necessary
condition, not a sufficient one: mRNA remains a proxy for protein activity,
and this script does not pretend otherwise.

Usage:
  python bench/analysis/empirical_axis_feasibility.py \
      --catalog <catalog dir> --curator-cases ab.tsv --experimental-cases exp.tsv \
      [--names ~/gitroot/mp-biopath-pathways/db_id_to_name_mapping.txt]
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import re
from collections import Counter

# A bare HGNC-style symbol and nothing else. Anything carrying a modification
# prefix (p-, 4xPalmC-, K6PolyUb-), a residue range (ACIN1(1-1093)) or a bare
# numeric id fails this deliberately.
BARE_SYMBOL = re.compile(r"^[A-Z][A-Z0-9-]{1,14}$")


def load_names(path: str) -> dict[str, str]:
    names: dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                names[parts[0].strip()] = parts[1].strip()
    return names


# Kinds that denote ONE entity, which is what "does this have a transcript"
# turns on. A dissociation_sink is a readout handle the generator mints for a
# released subunit -- it is still that single gene product, so it is grouped
# with simple_entity here. simple_complex and reaction are not single entities.
SINGLE_ENTITY_KINDS = ("simple_entity", "dissociation_sink")


def load_readout_kinds(catalog: str) -> dict[tuple[str, str], str]:
    """(pathway stable id, readout stable id) -> node_kind.

    Resolved PER PATHWAY, not catalog-wide. 2,785 of 15,616 stable ids in the
    92-pathway catalog appear as `simple_entity` in one pathway and
    `dissociation_sink` in another, so a catalog-wide first-match lookup is
    order-dependent: it silently reported a different answer depending on the
    order the directories were globbed in, and moved this script's headline by
    121 triples between two runs over the same data.
    """
    kinds: dict[tuple[str, str], str] = {}
    for d in sorted(glob.glob(os.path.join(catalog, "R-HSA-*"))):
        pid = os.path.basename(d)
        mapping = os.path.join(d, "stid_to_uuid_mapping.csv")
        nodes = os.path.join(d, "nodes.csv")
        if not (os.path.exists(mapping) and os.path.exists(nodes)):
            continue
        with open(nodes) as fh:
            uuid_kind = {r["uuid"]: r.get("node_kind", "") for r in csv.DictReader(fh)}
        with open(mapping) as fh:
            for r in csv.DictReader(fh):
                sid = str(r.get("stable_id", "")).replace("R-HSA-", "")
                kind = uuid_kind.get(r.get("uuid", ""), "")
                if not (sid and kind):
                    continue
                key = (pid, sid)
                # Within one pathway, prefer the entity reading over the sink
                # handle: the sink is a copy of the same species.
                if key not in kinds or (kinds[key] == "dissociation_sink"
                                        and kind == "simple_entity"):
                    kinds[key] = kind
    return kinds


def load_pathway_ids(path: str) -> dict[str, str]:
    """pathway benchmark name -> R-HSA id, from the generator pathway list."""
    out: dict[str, str] = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pid = str(r.get("id", "")).strip()
            nm = str(r.get("pathway_name", "")).strip()
            if pid and nm:
                out[nm] = pid if pid.startswith("R-HSA-") else f"R-HSA-{pid}"
    return out


def read_cases(path: str) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def triples(rows) -> set:
    return {(r["pathway"], r["gene"], r["key_output"]) for r in rows}


def symbol(names: dict[str, str], sid: str) -> str:
    return names.get(sid, "").split("_[")[0]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--curator-cases", required=True)
    ap.add_argument("--experimental-cases", required=True)
    ap.add_argument("--pathway-list", required=True,
                    help="generator pathway list (id, pathway_name) so readout kinds "
                         "can be resolved in the pathway the case belongs to")
    ap.add_argument("--names",
                    default=os.path.expanduser(
                        "~/gitroot/mp-biopath-pathways/db_id_to_name_mapping.txt"))
    ap.add_argument("--dump-genes", help="write the perturbed-gene and readout symbol lists here")
    a = ap.parse_args()

    names = load_names(a.names)
    kinds = load_readout_kinds(a.catalog)
    pid_of = load_pathway_ids(a.pathway_list)

    def kind_of(row) -> str:
        pid = pid_of.get(row["pathway"])
        if pid is None:
            return "unknown-pathway"
        return kinds.get((pid, row["key_output"]), "unmapped")
    curator = read_cases(a.curator_cases)
    experimental = read_cases(a.experimental_cases)

    covered = {r["pathway"] for r in experimental}
    held = [r for r in curator if r["pathway"] not in covered]

    n_curator_paths = len({r["pathway"] for r in curator})
    print(f"curator axis            {len(curator):>7,} cases  "
          f"{n_curator_paths:>3} pathways  {len(triples(curator)):>6,} triples")
    print(f"experimental axis today {len(experimental):>7,} cases  "
          f"{len(covered):>3} pathways  {len(triples(experimental)):>6,} triples   "
          f"(all inside the tuning ten -> no held-out empirical test)")
    print(f"pathways with curator truth and NO experimental truth: {len({r['pathway'] for r in held})}")
    print()
    print("=== funnel: what a transcript-based derived axis could score, held-out only ===")
    print(f"{'held-out curator cases':<46}{len(held):>7,} cases  {len(triples(held)):>6,} triples")

    by_kind = Counter(kind_of(r) for r in held)
    for kind, n in by_kind.most_common():
        print(f"    readout kind {kind:<28}{n:>7,} cases  ({100 * n / len(held):>4.1f}%)")

    single = [r for r in held if kind_of(r) in SINGLE_ENTITY_KINDS]
    print(f"{'single-entity readouts (entity or sink)':<46}{len(single):>7,} cases  {len(triples(single)):>6,} triples")

    clean = [r for r in single if BARE_SYMBOL.match(symbol(names, r["key_output"]))]
    print(f"{'  ... and a bare gene symbol (measurable)':<46}{len(clean):>7,} cases  {len(triples(clean)):>6,} triples")
    print()
    readouts = sorted({symbol(names, r["key_output"]) for r in clean})
    genes = sorted({r["gene"] for r in clean})
    paths = sorted({r["pathway"] for r in clean})
    print(f"ADDRESSABLE: {len(clean):,} cases / {len(triples(clean)):,} triples / "
          f"{len(readouts)} readout genes / {len(genes)} perturbed genes / {len(paths)} pathways")
    print("Every one of those pathways is outside the tuning ten, so this would be the")
    print("first held-out empirical test available to either tool.")
    print()
    print("REQUIRED EXTERNAL COVERAGE, the gating question:")
    print(f"  knockdown and over-expression signatures for {len(genes)} perturbed genes")
    print(f"  measured transcript levels for {len(readouts)} readout genes")
    print("  Caveat that no dataset removes: mRNA is a proxy for protein activity.")

    if a.dump_genes:
        with open(a.dump_genes, "w") as fh:
            fh.write("kind\tsymbol\n")
            for g in genes:
                fh.write(f"perturbed\t{g}\n")
            for r in readouts:
                fh.write(f"readout\t{r}\n")
        print(f"\nwrote gene lists to {a.dump_genes}")


if __name__ == "__main__":
    main()
