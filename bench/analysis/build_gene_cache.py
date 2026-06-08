#!/usr/bin/env python3
"""Pre-populate /tmp/gene_to_stids.json from Reactome's Neo4j.

The benchmark resolves gene names to Reactome PhysicalEntity stIds via
ReferenceEntity.geneName. We cache that mapping here so analysis scripts
don't need a live Neo4j connection.

Run once before any analysis. Re-run if new genes appear in the case dumps.
"""
import csv, json, sys
from pathlib import Path

DEFAULT_DUMP = "/tmp/ds_best_cases.tsv"
DEFAULT_CACHE = "/tmp/gene_to_stids.json"


def main(dump_path=DEFAULT_DUMP, cache_path=DEFAULT_CACHE):
    if not Path(dump_path).exists():
        print(f"ERROR: case dump not found at {dump_path}. "
              f"Run the benchmark with DS_DUMP_CASES=... first.", file=sys.stderr)
        sys.exit(2)
    genes = sorted({r["gene"] for r in csv.DictReader(
        open(dump_path), delimiter="\t")})
    print(f"{len(genes)} unique genes to resolve", file=sys.stderr)

    from py2neo import Graph
    g = Graph("bolt://localhost:7687", auth=("neo4j", "reactome"))
    rows = g.run(
        "UNWIND $names AS gn "
        "MATCH (re:ReferenceEntity)<-[:referenceEntity]-(pe:PhysicalEntity) "
        "WHERE gn IN re.geneName "
        "RETURN gn AS gene, COLLECT(DISTINCT pe.stId) AS stids",
        names=genes).data()
    gene_to_stids = {r["gene"]: list(r["stids"]) for r in rows}
    with open(cache_path, "w") as f:
        json.dump(gene_to_stids, f)
    print(f"Wrote {cache_path} with {len(gene_to_stids)} entries", file=sys.stderr)


if __name__ == "__main__":
    main(*sys.argv[1:])
