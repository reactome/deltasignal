"""Shared utilities for failure-analysis scripts.

Loads pathway networks from the upstream logic-network-generator catalog,
resolves gene → stids via a cached Neo4j lookup, and provides BFS / shortest-
path primitives on the logic network DAG.
"""
import csv, json, os
from collections import defaultdict, deque
from pathlib import Path

CATALOG = Path(os.environ.get(
    "PATHWAY_CATALOG",
    str(Path.home() / "gitroot" / "logic-network-generator" / "output")))

# Pathways used in the 9-pathway experimental benchmark (name → R-HSA id).
EXP_PATHWAYS = {
    "Transcriptional_Regulation_by_TP53": "3700989",
    "Mitotic_G1-G1_S_phases": "453279",
    "Signaling_by_ERBB2": "1227986",
    "Cell_Cycle_Checkpoints": "69620",
    "Signaling_by_WNT": "195721",
    "PIP3_activates_AKT_signaling": "1257604",
    "Mitotic_Prophase": "68875",
    "HDR_through_Homologous_Recombination_HRR_or_Single_Strand_Annealing_SSA_": "5693567",
    "S_Phase": "69242",
}

# Joint groups used by benchmark_selective_joint.py.
JOINT_GROUPS = {
    "Mitotic_G1-G1_S_phases":
        ["Mitotic_G1-G1_S_phases_R-HSA-453279", "S_Phase_R-HSA-69242"],
    "S_Phase":
        ["Mitotic_G1-G1_S_phases_R-HSA-453279", "S_Phase_R-HSA-69242"],
}

CLS_NAME = {"0": "DOWN", "1": "NORM", "2": "UP"}


def pw_dir(pname):
    """Pathway folder name in the catalog (handles trailing-underscore names)."""
    pid = EXP_PATHWAYS[pname]
    sep = "" if pname.endswith("_") else "_"
    return f"{pname}{sep}R-HSA-{pid}"


def load_one(dirname):
    """Load a single pathway. Returns (incoming, outgoing, stids, rev_stids, proxies)."""
    rows = list(csv.reader(open(CATALOG / dirname / "logic_network.csv")))[1:]
    incoming = defaultdict(list)
    outgoing = defaultdict(list)
    for r in rows:
        src, tgt, pn, ao, et, _ = r
        incoming[tgt].append((src, pn, ao, et))
        outgoing[src].append((tgt, pn, ao, et))
    stids, rev = {}, defaultdict(list)
    for r in csv.DictReader(open(CATALOG / dirname / "stid_to_uuid_mapping.csv")):
        stids[r["uuid"]] = r["stable_id"]
        rev[r["stable_id"]].append(r["uuid"])
    proxies = defaultdict(list)
    pf = CATALOG / dirname / "entity_reaction_proxy_mapping.csv"
    if pf.exists():
        for r in csv.DictReader(open(pf)):
            proxies[r["entity_stable_id"]].append(r["proxy_uuid"])
    return incoming, outgoing, stids, dict(rev), dict(proxies)


def load_joint(dirnames):
    """Union the logic networks of multiple pathways. UUIDs are unique per pathway,
    so the union is well-formed without renaming."""
    inc, out = defaultdict(list), defaultdict(list)
    stids, rev, proxies = {}, defaultdict(list), defaultdict(list)
    for d in dirnames:
        rows = list(csv.reader(open(CATALOG / d / "logic_network.csv")))[1:]
        for r in rows:
            src, tgt, pn, ao, et, _ = r
            inc[tgt].append((src, pn, ao, et))
            out[src].append((tgt, pn, ao, et))
        for r in csv.DictReader(open(CATALOG / d / "stid_to_uuid_mapping.csv")):
            stids[r["uuid"]] = r["stable_id"]
            rev[r["stable_id"]].append(r["uuid"])
        pf = CATALOG / d / "entity_reaction_proxy_mapping.csv"
        if pf.exists():
            for r in csv.DictReader(open(pf)):
                proxies[r["entity_stable_id"]].append(r["proxy_uuid"])
    return inc, out, stids, dict(rev), dict(proxies)


class Networks:
    """Lazy-load and cache per-pathway / per-joint networks."""

    def __init__(self):
        self._single = {}
        self._joint = {}

    def get(self, pname):
        """Return the network the benchmark would have used for this pathway:
        joint if the pathway participates in a JOINT_GROUPS entry, else single."""
        if pname in JOINT_GROUPS:
            if pname not in self._joint:
                self._joint[pname] = load_joint(JOINT_GROUPS[pname])
            return self._joint[pname]
        if pname not in self._single:
            self._single[pname] = load_one(pw_dir(pname))
        return self._single[pname]

    def get_single(self, pname):
        """Always return the SINGLE pathway network (ignores joint)."""
        if pname not in self._single:
            self._single[pname] = load_one(pw_dir(pname))
        return self._single[pname]


def load_gene_stids(cache_path="/tmp/gene_to_stids.json"):
    """Returns dict gene_name → set of stable_ids (from Reactome ReferenceEntity)."""
    if not os.path.exists(cache_path):
        raise SystemExit(
            f"No gene→stids cache at {cache_path}. "
            f"Run `python bench/analysis/build_gene_cache.py` first.")
    return {k: set(v) for k, v in json.load(open(cache_path)).items()}


def bfs_upstream(incoming, sinks, max_depth=30):
    """Return all nodes upstream (backward) of any sink, within max_depth hops."""
    seen = set(sinks)
    q = deque((s, 0) for s in sinks)
    while q:
        u, d = q.popleft()
        if d >= max_depth:
            continue
        for src, _, _, _ in incoming.get(u, []):
            if src not in seen:
                seen.add(src)
                q.append((src, d + 1))
    return seen


def shortest_path_back(incoming, sources, sinks, max_depth=30):
    """BFS from any sink backward until reaching any source.
    Returns a list of edge tuples (src, tgt, pn, ao, et), src→...→sink.
    None if no path exists."""
    parents = {s: None for s in sinks}
    q = deque(sinks)
    found = None
    while q:
        u = q.popleft()
        for src, pn, ao, et in incoming.get(u, []):
            if src not in parents:
                parents[src] = (u, pn, ao, et)
                if src in sources:
                    found = src
                    q.clear()
                    break
                q.append(src)
    if found is None:
        return None
    edges = []
    cur = found
    while parents[cur] is not None:
        nxt, pn, ao, et = parents[cur]
        edges.append((cur, nxt, pn, ao, et))
        cur = nxt
    return edges


def resolve_gene_uuids(gene, gene_stids, rev_stids):
    """Gene name → list of UUIDs in this pathway network."""
    return [u for s in gene_stids.get(gene, set()) for u in rev_stids.get(s, [])]


def resolve_readout_uuids(key_output, rev_stids, proxies):
    """key_output (numeric reactome id) → list of UUIDs in this network.
    Falls back to proxy mapping if the readout is decomposed into producing reactions."""
    ko_stid = f"R-HSA-{key_output}"
    return rev_stids.get(ko_stid, []) or proxies.get(ko_stid, [])


def read_case_dump(path):
    """Read a per-case TSV dump (produced by DS_DUMP_CASES in the benchmark)."""
    return list(csv.DictReader(open(path), delimiter="\t"))
