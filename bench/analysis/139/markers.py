#!/usr/bin/env python3
"""specs/139: each CPPA drug's target-activity markers, from Reactome Release 97 and the drug targets (values-blind).

A marker of drug d is a CPPA phospho antibody with a site Reactome curates as made by a reaction (the phosphorylated
protein in an output, to depth 3, and in no input) catalysed by one of d's primary targets: the catalyst activity's
activeUnit if any, else its physical entity, counted only when that unit is not a complex. Primary targets:
kinobeads targets within 10-fold of the drug's lowest Kd_app, else ChEMBL mechanism targets, else an anonymised
compound's class genes; ChEMBL binding alone gives none. Every query is bounded on the server
(apoc.cypher.runTimeboxed), so an interrupted client cannot leave it running.

Usage: markers.py READOUTS_TSV DRUG_TARGETS_TSV ACC_GENE_TSV > markers.tsv   (the output stays outside the repository)
"""
import collections
import csv
import subprocess
import sys

TIMEBOX_MS = 240000
END = "__end_of_rows__"
KD_FOLD = 10
SHELL = ["docker", "exec", "-i", "reactome-neo4j", "sh", "-c",
         'cypher-shell -u neo4j -p "${NEO4J_AUTH#neo4j/}" --format plain']


def cypher(inner, columns):
    """Rows of a read query run inside apoc.cypher.runTimeboxed; `inner` must not contain double quotes. The timebox
    returns the rows streamed so far when it stops a query, and swallows an inner error as zero rows, so a sentinel row
    is appended (UNION ALL streams it last) and its absence is an error, not a short answer."""
    end = ", ".join(f"'{END}' AS {c}" for c in columns)
    q = f'CALL apoc.cypher.runTimeboxed("{inner} UNION ALL RETURN {end}", {{}}, {TIMEBOX_MS}) YIELD value RETURN ' + \
        ", ".join(f"value.{c} AS {c}" for c in columns) + ";"
    out = subprocess.run(SHELL, input=q, capture_output=True, text=True, timeout=TIMEBOX_MS // 1000 + 60, check=True)
    rows = [[f.strip().strip('"') for f in line.split(", ")] for line in out.stdout.splitlines()[1:] if line]
    if not rows or rows[-1] != [END] * len(columns):
        raise SystemExit(f"Neo4j query incomplete (timed out or failed inside the timebox): {inner[:120]}...")
    return rows[:-1]


def quoted(values):
    return ",".join(f"'{v}'" for v in values)


def phospho_sites(readouts_tsv):
    """{antibody: ['GENE:coord', ...]} for phospho readouts."""
    out = {}
    for line in open(readouts_tsv):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) > 2 and f[1] == "phospho":
            out[f[0]] = f[2].split(",")
    return out


def site_catalysts(sites, gene_accs):
    """{'GENE:coord': {accession of a non-complex catalyst unit of a reaction making that phospho site}}."""
    accs = sorted({a for s in sites for a in gene_accs.get(s.split(":")[0], ())})
    pairs = ",".join(f"['{s.split(':')[0]}',{s.split(':')[1]}]" for s in sorted(sites))
    ewas = cypher(
        f"MATCH (rgp:ReferenceGeneProduct) WHERE rgp.identifier IN [{quoted(accs)}] "
        "MATCH (rgp)<-[:referenceEntity]-(e:EntityWithAccessionedSequence)-[:hasModifiedResidue]->(mr)-[:psiMod]->(pm) "
        "WHERE pm.name[0] CONTAINS 'phospho' "
        "WITH e, [g IN rgp.geneName | g + ':' + toString(mr.coordinate)] AS keys "
        f"UNWIND [{pairs}] AS gs WITH e, gs, keys WHERE (gs[0] + ':' + toString(gs[1])) IN keys "
        "RETURN DISTINCT gs[0] + ':' + toString(gs[1]) AS site, e.stId AS ewas", ["site", "ewas"])
    rx = cypher(
        "UNWIND [" + ",".join(f"['{s}','{e}']" for s, e in ewas) + "] AS p "
        "MATCH (e:PhysicalEntity {stId: p[1]}) "
        "MATCH (e)<-[:hasComponent|hasMember|hasCandidate*0..3]-(:PhysicalEntity)<-[:output]-(r:ReactionLikeEvent) "
        "WHERE NOT EXISTS { MATCH (r)-[:input]->(:PhysicalEntity)-[:hasComponent|hasMember|hasCandidate*0..3]->(e) } "
        "RETURN DISTINCT p[0] AS site, r.stId AS reaction", ["site", "reaction"])
    units = cypher(
        f"UNWIND [{quoted(sorted({r for _, r in rx}))}] AS rid "
        "MATCH (r:ReactionLikeEvent {stId: rid})-[:catalystActivity]->(ca:CatalystActivity)-[:physicalEntity]->(c) "
        "OPTIONAL MATCH (ca)-[:activeUnit]->(au:PhysicalEntity) "
        "WITH rid, CASE WHEN au IS NULL THEN c ELSE au END AS unit WHERE NOT unit:Complex "
        "MATCH (unit)-[:hasMember|hasCandidate*0..3]->(:EntityWithAccessionedSequence)-[:referenceEntity]->(g) "
        "RETURN DISTINCT rid AS reaction, g.identifier AS acc", ["reaction", "acc"])
    by_rx = collections.defaultdict(set)
    for r, a in units:
        by_rx[r].add(a)
    out = collections.defaultdict(set)
    for s, r in rx:
        out[s] |= by_rx[r]
    return out


def primary_targets(rows, gene_accs):
    """A drug's primary target accessions from its drug_targets.tsv rows (see the module docstring)."""
    src = {r["source"] for r in rows}
    if "kinobeads" in src:
        kd = [(float(r["kd_app_nM"]), r["target"]) for r in rows if r["source"] == "kinobeads"]
        best = min(k for k, _ in kd)
        return {t for k, t in kd if k <= KD_FOLD * best}
    out = set()
    for r in rows:
        if r["source"] == "class":
            out |= gene_accs.get(r["target"], set())
        elif r["source"] == "mechanism":
            out.add(r["target"])
    return out


def markers(antibody_sites, catalysts, targets):
    return sorted(a for a, ss in antibody_sites.items() if any(targets & catalysts.get(s, set()) for s in ss))


def main():
    readouts_tsv, targets_tsv, acc_gene_tsv = sys.argv[1:]
    gene_accs = collections.defaultdict(set)
    for r in csv.DictReader(open(acc_gene_tsv), delimiter="\t"):
        gene_accs[r["gene"]].add(r["acc"])
    ab = phospho_sites(readouts_tsv)
    cat = site_catalysts({s for ss in ab.values() for s in ss}, gene_accs)
    by_drug = collections.defaultdict(list)
    for r in csv.DictReader(open(targets_tsv), delimiter="\t"):
        by_drug[r["drug"]].append(r)
    print("drug\trole\tprimary_targets\tmarkers")
    for d in sorted(by_drug):
        t = primary_targets(by_drug[d], gene_accs)
        print(f"{d}\t{by_drug[d][0]['role']}\t{','.join(sorted(t))}\t{','.join(markers(ab, cat, t))}")


if __name__ == "__main__":
    main()
