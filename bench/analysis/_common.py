"""Shared utilities for failure-analysis scripts.

Loads pathway networks from the upstream logic-network-generator catalog,
resolves gene → stids via a cached Neo4j lookup, and provides BFS / shortest-
path primitives on the logic network DAG.
"""
import csv
import json
import os
import re
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
        src, tgt, pn, ao, et, *_ = r   # 6 or 7 columns (edge_reaction_id added 2026-09)
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
            src, tgt, pn, ao, et, *_ = r   # 6 or 7 columns (edge_reaction_id added 2026-09)
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

def pathway_dir_index(root: Path) -> dict[str, Path]:
    """Map Reactome numeric id -> catalog directory.

    Keyed by ID, never by pathway name. Names are data: they gain commas, lose
    trailing underscores and get recurated, and three of them differed between
    the curator files and the catalog on Release97 -- which silently moved 438
    of 1,484 cases into a bogus bucket before it was noticed.

    Accepts both the historical `Some_Pathway_Name_R-HSA-12345` layout and a
    bare `R-HSA-12345`, so it keeps working across the naming change.

    Matches any species prefix, not just HSA. The generator names a directory by
    whatever stable id it is given, so an `R-MMU-` directory is possible; an
    HSA-only pattern would skip it SILENTLY, which is the same failure mode this
    function exists to remove. Keyed by the numeric part, because the case dumps
    and pathway lists carry numeric ids — so a cross-species catalog could
    collide, and this asserts rather than picking one at random.
    """
    out: dict[str, Path] = {}
    for d in root.iterdir():
        if not d.is_dir():
            continue
        m = re.search(r"R-[A-Z]{3}-(\d+)$", d.name)
        if not m:
            continue
        pid = m.group(1)
        if pid in out:
            raise ValueError(
                f"Two directories share pathway id {pid}: {out[pid].name} and "
                f"{d.name}. Ids are keyed numerically here; a mixed-species "
                "catalog needs the full stable id as the key."
            )
        out[pid] = d
    return out


def name_to_id_map(pathway_list: Path) -> dict[str, str]:
    """Curator pathway NAME -> numeric id, tolerating the known spelling drift.

    Only for reading legacy case dumps, which record a name and not an id.
    Anything writing new data should carry the id.
    """
    raw: dict[str, str] = {}
    with pathway_list.open(newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            raw[r["pathway_name"]] = r["pathway_id"]
    norm = {k.replace(",", "").rstrip("_"): v for k, v in raw.items()}

    def lookup(name: str) -> str | None:
        return raw.get(name) or norm.get(name.replace(",", "").rstrip("_"))

    return type("M", (), {"get": staticmethod(lambda n, d=None: lookup(n) or d)})()


def mcnemar_exact(fixed: int, broke: int) -> float:
    """Two-sided exact McNemar p-value for a paired arm comparison.

    `fixed` and `broke` are the DISCORDANT pairs: cases the arm got right and
    the baseline got wrong, and vice versa. Concordant pairs carry no
    information about which arm is better and are excluded, which is the whole
    point of the test.

    This exists because net case deltas were being reported as though their
    sign were meaningful. They are not, on their own: an AND-mode arm read
    "-2 of 849" on the experimental axis, which is 0 fixed and 2 broke and
    p = 0.50 -- pure noise -- while "-16 of 18,808" on the curator axis is 11
    fixed and 30 broke and p = 0.0043, a real signal. The two look comparable
    as percentages and are not comparable at all.

    Report the p-value beside every net figure, or the reader cannot tell a
    result from a coin flip.
    """
    from math import comb

    n = fixed + broke
    if n == 0:
        return 1.0
    k = min(fixed, broke)
    tail = sum(comb(n, i) for i in range(k + 1)) / 2 ** n
    return min(1.0, 2 * tail)
