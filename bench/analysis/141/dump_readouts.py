#!/usr/bin/env python3
"""specs/141 FR-003/FR-004, bench side: dump what the BENCH resolves for every
benchmark readout, via the very functions the benchmarks run, so the engine's
/api/readouts port can be held equal to them.

Two sections:
- curator: every (pathway, key_output) of the curator benchmark, resolved by
  benchmark_vs_mpbiopath.resolve_key_output (the ladder run_pathway calls),
  with the candidate stable ids, the basis, the curated-entity groups, and
  the Neo4j member list the set branch consulted (fed to the engine, which
  has no Neo4j by design).
- korkut: every total readout of the korkut bench on every pathway it is
  resolved in, via korkut_bench.readout_nodes and abundance_copies
  (abundance_a3, the adopted rule): the stable ids, the resolved occurrence
  list, and the kept abundance copies.

Run with DS_CATALOG_ROOT pointing at the catalog under test (the frozen build
for the registered parity run):
  DS_CATALOG_ROOT=~/deltasignal-catalogs/current \
    python bench/analysis/141/dump_readouts.py bench/analysis/141/bench_readouts.json
"""
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # the bench dir

import benchmark_vs_mpbiopath as bm  # noqa: E402
import korkut_bench as kb  # noqa: E402


def pathway_list():
    with open(bm.MPBIO_ROOT / "pathway_list.tsv") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx_id, idx_name = header.index("pathway_id"), header.index("pathway_name")
        pathways = [(c[idx_id], c[idx_name]) for c in
                    (line.rstrip("\n").split("\t") for line in f)
                    if len(c) > max(idx_id, idx_name)]
    return bm.apply_catalog_ids(pathways,
                                Path(__file__).resolve().parents[2] / "catalog_pathways.tsv")


def dump_curator():
    out = []
    for pid, pname in pathway_list():
        pathway_dir = bm.find_pathway_dir(pid)
        if pathway_dir is None:
            continue
        curator = bm.load_curator(pname)
        if curator is None:
            continue
        _header, rows = curator
        stid_to_uuids = bm.load_stid_to_uuids(pathway_dir)
        proxies = bm.load_entity_reaction_proxies(pathway_dir)
        output_entities = bm.load_output_entities(pathway_dir)
        own_uuids = bm.context_own_uuids(pathway_dir)
        kos = []
        seen = set()
        for r in rows:
            ko = str(r["key_output"])
            if not ko or ko in seen:
                continue
            seen.add(ko)
            kos.append(ko)
        dbid_to_stid = bm.neo4j_dbid_to_stid(set(kos))
        readouts = []
        for ko in kos:
            candidate_sids = []
            resolved = dbid_to_stid.get(ko)
            if resolved:
                candidate_sids.append(resolved)
            if f"R-HSA-{ko}" not in candidate_sids:
                candidate_sids.append(f"R-HSA-{ko}")
            _uuids, groups, basis, members = bm.resolve_key_output(
                candidate_sids, stid_to_uuids, proxies, output_entities,
                own_uuids=own_uuids)
            readouts.append({"ko": ko, "stable_ids": candidate_sids,
                             "basis": basis, "groups": groups, "members": members})
        out.append({"pathway_id": pid, "pathway_name": pname,
                    "dir": pathway_dir.name, "readouts": readouts})
        print(f"curator {pname}: {len(readouts)} readouts", file=sys.stderr, flush=True)
    return out


def dump_korkut():
    if kb.TOTAL_READOUT not in kb.ABUNDANCE_RULES:
        raise SystemExit(f"DS_TOTAL_READOUT={kb.TOTAL_READOUT}: the gate runs the adopted abundance rule")
    readouts = [r for r in kb.read_tsv(kb.MAPS / "readouts.tsv")
                if r["kind"] == "total"]
    genes = sorted({g for r in readouts for g in r["sites"].split(",")})
    gene_to_stids = bm.neo4j_gene_to_stids(genes)
    r_stids = kb.readout_stids(readouts, kb.total_gene_stids(readouts, gene_to_stids), {})
    out = []
    build = bm.CATALOG_ROOT.resolve()
    for pdir in sorted(p for p in build.iterdir() if (p / "stid_to_uuid_mapping.csv").exists()):
        stid_to_uuids = bm.load_stid_to_uuids(pdir)
        ro_nodes = kb.readout_nodes(r_stids, stid_to_uuids, bm.context_own_uuids(pdir))
        if not ro_nodes:
            continue
        with bm.open_network(pdir) as f:
            net_rows = list(csv.DictReader(f))
        entries = []
        for name in sorted(ro_nodes):
            nodes = ro_nodes[name]
            kept = kb.abundance_copies(nodes, net_rows, kb.TOTAL_READOUT)
            entries.append({"readout": name, "stable_ids": r_stids[name],
                            "occurrences": nodes, "kept": kept})
        out.append({"dir": pdir.name, "readouts": entries})
        print(f"korkut {pdir.name}: {len(entries)} total readouts", file=sys.stderr, flush=True)
    return out


def main():
    out_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / "bench_readouts.json"
    print(f"catalog: {bm.CATALOG_ROOT}", flush=True)
    data = {"catalog": bm.CATALOG_ROOT.resolve().name,
            "abundance_rule": kb.TOTAL_READOUT,
            "curator": dump_curator(),
            "korkut": dump_korkut()}
    n_cur = sum(len(p["readouts"]) for p in data["curator"])
    n_kor = sum(len(p["readouts"]) for p in data["korkut"])
    out_path.write_text(json.dumps(data, indent=1))
    print(f"wrote {out_path}: {n_cur} curator readout resolutions over {len(data['curator'])} pathways, "
          f"{n_kor} korkut total readouts over {len(data['korkut'])} pathways", flush=True)


if __name__ == "__main__":
    main()
