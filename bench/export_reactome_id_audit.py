#!/usr/bin/env python3
"""Export benchmark-relevant Reactome identifiers from a pinned Neo4j release."""

from __future__ import annotations

import argparse
import base64
import csv
import json
from pathlib import Path
from urllib.request import Request, urlopen

from benchmark_mpbiopath_cases import load_cases, load_gene_dbids


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--supplementary-workbook", type=Path, required=True)
    parser.add_argument("--id-map", type=Path, required=True)
    parser.add_argument(
        "--neo4j-http",
        default="http://127.0.0.1:7474/db/graph.db/tx/commit",
    )
    parser.add_argument("--neo4j-user", default="neo4j")
    parser.add_argument("--neo4j-password", default="")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=500)
    return parser.parse_args()


def query_batch(
    endpoint: str,
    user: str,
    password: str,
    dbids: list[int],
) -> list[dict[str, str]]:
    query = """
    MATCH (n)
    WHERE n.dbId IN $dbids
    RETURN toString(n.dbId) AS dbid,
           coalesce(n.stId, '') AS stable_id,
           coalesce(n.schemaClass, head(labels(n)), '') AS schema_class,
           coalesce(n.displayName, '') AS display_name
    ORDER BY n.dbId
    """
    payload = json.dumps(
        {"statements": [{"statement": query, "parameters": {"dbids": dbids}}]}
    ).encode()
    headers = {"Content-Type": "application/json"}
    if user:
        token = base64.b64encode(f"{user}:{password}".encode()).decode()
        headers["Authorization"] = f"Basic {token}"
    request = Request(endpoint, data=payload, headers=headers)
    with urlopen(request, timeout=120) as response:
        body = json.loads(response.read())
    errors = body.get("errors", [])
    if errors:
        raise RuntimeError(f"Neo4j query failed: {errors}")
    result = body["results"][0]
    columns = result["columns"]
    return [dict(zip(columns, item["row"])) for item in result["data"]]


def main() -> None:
    args = parse_args()
    cases = load_cases(args.supplementary_workbook, "experimental")
    genes = {case.gene for case in cases}
    gene_dbids = load_gene_dbids(args.id_map, genes)
    requested = {int(case.key_output_dbid) for case in cases}
    requested.update(
        int(dbid)
        for identifiers in gene_dbids.values()
        for dbid in identifiers
        if dbid.isdigit()
    )

    rows: list[dict[str, str]] = []
    ordered = sorted(requested)
    for start in range(0, len(ordered), args.batch_size):
        rows.extend(
            query_batch(
                args.neo4j_http,
                args.neo4j_user,
                args.neo4j_password,
                ordered[start : start + args.batch_size],
            )
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("dbid", "stable_id", "schema_class", "display_name"),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)
    present = {int(row["dbid"]) for row in rows}
    print(
        json.dumps(
            {
                "requested_dbids": len(requested),
                "present_in_release": len(present),
                "absent_from_release": len(requested - present),
                "output": str(args.output),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
