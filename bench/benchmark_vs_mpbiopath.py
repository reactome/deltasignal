#!/usr/bin/env python3
"""DeltaSignal vs naive-Boolean benchmark against the MP-BioPath curator set.

Runs the same protocol the logic-network-generator uses in
`bin/validate-against-mpbiopath.py`, but instead of the 3-valued Boolean
propagator we call deltasignal's /api/solve and compare the resulting
continuous node activities back to the curator's discrete {0,1,2} state.

Headline numbers we're aiming to beat:
  - MP-BioPath original tool on v86 networks:  ~75% vs experimental
  - Naive Boolean on regenerated v96 networks: 70.55% vs curator (12,895 cases)
  - Curator vs experimental:                   ~81% (their human ceiling)

Discretization for deltasignal:
  Pinned perturbation: direction=0 → activity=0,  direction=2 → activity=80
  (DS_PERTURB_UI_DOWN / DS_PERTURB_UI_UP override these, for the specs/021
  dose-response runs only).
  Read key-output's solved UI value (max over its UUIDs).
  Classify:
    < 0.5  → DOWN (0)
    [0.5, 2.0) → NORMAL (1)   (baseline UI = 1)
    >= 2.0 → UP (2)
"""

import argparse
import csv
import json
import os
import re
import sys
import time
from collections import defaultdict, Counter
from collections.abc import Mapping
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

# Locations of the upstream data, overridable via env so the harness isn't
# pinned to one developer's checkout:
#   LOGIC_NETWORK_GENERATOR — the logic-network-generator repo (with output/)
#   MPBIO_PATHWAYS          — the mp-biopath-pathways ground-truth repo
GENERATOR_ROOT = Path(os.environ.get(
    "LOGIC_NETWORK_GENERATOR",
    str(Path.home() / "gitroot" / "logic-network-generator")))
MPBIO_ROOT = Path(os.environ.get(
    "MPBIO_PATHWAYS",
    str(Path.home() / "gitroot" / "mp-biopath-pathways")))
# The generator's own output/ by default, but overridable: a validation run
# needs to point at a specific catalog build, and symlinking around a
# hard-coded path is how you end up scoring the wrong networks.
CATALOG_ROOT = Path(os.environ.get("DS_CATALOG_ROOT", str(GENERATOR_ROOT / "output")))
CURATOR_DIR = MPBIO_ROOT / "reactome_curator_predictions"
EXPERIMENTAL_DIR = MPBIO_ROOT / "experimental_results"

# DeltaSignal HTTP API
DS_BASE = os.environ.get("DELTASIGNAL_BASE", "http://127.0.0.1:8080")

# Neo4j (read-only; uses the same image the generator uses)
NEO4J_URL = os.environ.get("NEO4J_URL", "bolt://localhost:7687")
NEO4J_USER = os.environ.get("NEO4J_USER", "neo4j")
NEO4J_PASSWORD = os.environ.get("NEO4J_PASSWORD", "")

DOWN, NORMAL, UP = 0, 1, 2

# Discretization thresholds on the UI 0-100 scale where 1 = baseline.
# Conservative thresholds: an unperturbed node holds at ~1.0 and is normal;
# a knockout that propagates correctly should land well below 0.5;
# an upregulation should clear 2.0.
DOWN_CUTOFF = float(os.environ.get("DS_DOWN_CUTOFF", 0.85))  # UI < cutoff → DOWN
UP_CUTOFF = float(os.environ.get("DS_UP_CUTOFF", 1.15))      # UI >= cutoff → UP

# Perturbation strength. The defaults are what every benchmark number in this
# repo was measured at. They are overridable only so a dose-response study can
# rerun the IDENTICAL pipeline at other strengths (specs/021); the curator and
# experimental ground truths are defined at full strength, so accuracy figures
# at any other value are not comparable to the published ones.
PERTURB_UI_DOWN = float(os.environ.get("DS_PERTURB_UI_DOWN", 0.0))   # direction=0: knockout
PERTURB_UI_UP = float(os.environ.get("DS_PERTURB_UI_UP", 80.0))      # direction=2: strong upregulation
if not (0.0 <= PERTURB_UI_DOWN < 1.0 < PERTURB_UI_UP <= 100.0):
    raise SystemExit(
        f"DS_PERTURB_UI_DOWN={PERTURB_UI_DOWN} and DS_PERTURB_UI_UP={PERTURB_UI_UP} must "
        "satisfy 0 <= down < 1 (baseline) < up <= 100; a 'knockdown' above baseline or an "
        "'overexpression' below it would silently invert the perturbation.")
PIN_CONFIDENCE = 1.0

# How to collapse a key-output's multiple UUID activities into one prediction.
KO_AGG = os.environ.get("DS_KO_AGG", "max")

# Diagnostic: comma-separated logic_network edge_types to drop before solving
# (e.g. "assembly,dissociation"). Empty = keep all edges.
SKIP_EDGE_TYPES = {t.strip() for t in os.environ.get("DS_SKIP_EDGE_TYPES", "").split(",") if t.strip()}
# specs/025 (Adam's rule, 2026-09-25). A regulator that names a SPECIFIC member
# of a set -- CREBBP:NS1, not the CREBBP/EP300 set -- regulates only the
# variant reactions built from that member. The generator attaches a Reactome
# reaction's regulators to every variant, so knocking CREBBP down removed
# CREBBP:NS1 from the EP300 variant too and de-repressed it 10x.
#   off (default) | member_only
SIBLING_REGULATORS = os.environ.get("DS_SIBLING_REGULATORS", "off")
if SIBLING_REGULATORS not in ("off", "member_only"):
    raise SystemExit(f"DS_SIBLING_REGULATORS={SIBLING_REGULATORS!r} must be 'off' or 'member_only'")
# Diagnostic: drop inhibitor edges whose source CONTAINS one of the same
# reaction's activators. Such an inhibitor is not an independent variable --
# it is a partition of the substrate pool, and it rises because the substrate
# rose. With divide-form inhibition (H = baseline/x) the two cancel, so an
# elevated substrate produces exactly baseline flux.
#
# Traced from a real failure: ATM over-expression in Cell_Cycle_Checkpoints
# drives CDKN1A mRNA to 76.83x, and "PCBP4 modulates CDKN1A translation" has
# that mRNA as its activator (76.83x) and PCBP4:CDKN1A mRNA -- a complex
# CONTAINING the mRNA -- as its inhibitor, also 76.83x. 76.83 * (0.01/0.7683)
# = 0.9998, so the readout reads exactly baseline and 10 steps of correct
# propagation are annulled.
#
# 150 reactions catalog-wide, 22.7% of those carrying both an activator and an
# inhibitor, across 37 pathways.
SKIP_SELF_CONTAINED_INHIBITORS = os.environ.get("DS_SKIP_SELF_INH", "0") == "1"
# Which of the nodes a gene resolves to get pinned. specs/023.
#   root  (the protocol of record's base -- Adam, 2026-09-25, following the
#         MP-BioPath publication): only ROOT inputs (no incoming edge) that are
#         or contain the gene. Everything downstream is computed. A gene with no
#         root form is not perturbed, and its cases are invalid, not scored.
#   entry: where the gene ENTERS the network -- resolved nodes no other resolved
#         node reaches, so a gene with no root form gets its first occurrence.
#   all   (the protocol from 0bd4565, 2026-07-14, to 2026-09-25): every node
#         whose members include the gene, wherever it sits -- 89% of pins were
#         mid-pathway complexes, SET rather than computed. Kept only so older
#         results can be reproduced.
#   root_cycle (DEFAULT since 2026-09-25, specs/024): `root`, but a gene whose only forms sit inside a catalytic
#         cycle -- the enzyme is regenerated by its own reaction, so it is fed
#         only by its own downstream (or our derived dissociation/depletion
#         edges) -- has that recycled form pinned instead. MP-BioPath's
#         published networks had those loops cut by hand, leaving the free
#         enzyme as a ROOT (specs/024; ALKBH2, CBLB traced). Genes with a true
#         root are unchanged.
PIN_SCOPE = os.environ.get("DS_PIN_SCOPE", "root_cycle")   # default since specs/024 (Adam)
PIN_TALLY: Counter = Counter()
if PIN_SCOPE not in ("all", "entry", "root", "root_cycle"):
    raise SystemExit(f"DS_PIN_SCOPE={PIN_SCOPE!r} must be 'root', 'root_cycle', 'entry' or 'all'")
# Diagnostic: collapse duplicate ACTIVATOR edges from the same source into the
# same reaction. An entity that is both the catalyst and a substrate of one
# reaction currently contributes TWICE to the AND product, so its fold-change
# is squared -- a 5x input yields 25x.
#
# Traced from a real failure: TP53 knockout in Cell_Cycle_Checkpoints reaches
# "MDM2 ubiquitinates phosphorylated MDM4", whose activator edges are Ub (1x)
# plus MDM2:MDM4 as BOTH input and catalyst (5x each). 5 * 5 = 25, and the
# readout reads UP where the curator expects no change.
#
# 7,013 reactions catalog-wide -- 15.8% of all reactions with 2+ activator
# edges, across 78 pathways -- and the repeated role pair is catalyst+input in
# 7,042 of 7,045 cases.
#
# NOTE: a "dedup activators" arm was recorded as negative previously, at
# p = 0.25 on the 742-case set. That is not significant, and that set has
# reversed decisions twice. This re-measures on the wide curator set.
DEDUP_ACTIVATORS = os.environ.get("DS_DEDUP_ACTIVATORS", "0") == "1"

# Key-output remap for pathways whose 2019 curator ground truth references
# entities that were deleted/replaced by a later Reactome recuration (the old
# dbId no longer resolves, so the case is unscoreable). Opt-in via DS_KO_REMAP=1.
#
# Currently handles the RHO GTPase cycle: the 2019 model had two GENERIC pooled
# readouts — RhoGTPase:GTP (dbId 194890, active) and RhoGTPase:GDP (194900,
# inactive). v97 replaced these with PER-PARALOG forms (RHOA:GTP, RAC1:GTP, …).
# We remap the generic readout to the perturbed gene's OWN active/inactive form
# when it is itself a GTPase, else to the pooled union of all such forms in the
# network (matching the original pooled semantics for GEF/GAP/GDI perturbations).
KO_REMAP = os.environ.get("DS_KO_REMAP", "0") == "1"
GENERIC_GTPASE_READOUTS = {"194890": "GTP", "194900": "GDP"}


def classify(ui_value: float) -> int:
    if ui_value < DOWN_CUTOFF:
        return DOWN
    if ui_value >= UP_CUTOFF:
        return UP
    return NORMAL


def apply_catalog_ids(pathways: list, catalog_list: Path) -> list:
    """Take each pathway's id from bench/catalog_pathways.tsv, the list the
    catalog is built from, instead of MP-BioPath's pathway_list.tsv. They must
    agree, and did not: MP-BioPath lists Interleukin-2_family_signaling as
    447115, which is Interleukin-12 family signaling (specs/024). With the
    catalog corrected and this unchanged, the benchmark looked for 447115,
    found no network, and silently dropped all 260 IL-2 family cases."""
    if not catalog_list.exists():
        return pathways
    ids = {}
    for line in open(catalog_list):
        if line.startswith("#") or line.startswith("id\t") or not line.strip():
            continue
        sid, name = line.rstrip("\n").split("\t")
        ids[name] = sid.rsplit("-", 1)[-1]
    out = []
    for pid, name in pathways:
        new = ids.get(name, pid)
        if new != pid:
            print(f"  [id] {name}: {pid} -> {new} (from {catalog_list.name})", flush=True)
        out.append((new, name))
    return out


def find_pathway_dir(numeric_id: str):
    """The R-HSA-suffixed dir is preferred (see generator dedupe convention)."""
    for d in sorted(CATALOG_ROOT.iterdir()):
        if d.is_dir() and d.name.endswith(f"R-HSA-{numeric_id}"):
            return d
    return None


def load_stid_to_uuids(pathway_dir: Path):
    """stable_id → list of UUIDs.

    Variant-aware: when the generator emits ``nodes.csv`` (complexes bundled as
    set_variant nodes), a plain Reactome stId no longer appears verbatim as a
    node id. We index each node by BOTH its ``diagram_entity_id`` (the stId a
    diagram renders — parent complex for a variant) AND every stId in its
    ``member_leaves``. So a gene resolves to every variant node that contains
    it, and a complex/entity readout resolves to its node(s) — no string
    parsing. Falls back to the legacy stid_to_uuid_mapping.csv when nodes.csv is
    absent (older catalogs).
    """
    out = defaultdict(list)
    nodes_csv = pathway_dir / "nodes.csv"
    if nodes_csv.exists():
        seen = set()
        with open(nodes_csv) as f:
            for row in csv.DictReader(f):
                uuid = str(row["uuid"])
                keys = set()
                de = (row.get("diagram_entity_id") or "").strip()
                if de:
                    keys.add(de)
                for m in (row.get("member_leaves") or "").split("|"):
                    m = m.strip()
                    if m:
                        keys.add(m)
                for k in keys:
                    if (k, uuid) not in seen:
                        seen.add((k, uuid))
                        out[k].append(uuid)
        return out
    with open(pathway_dir / "stid_to_uuid_mapping.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            out[str(row["stable_id"])].append(str(row["uuid"]))
    return out


def load_entity_reaction_proxies(pathway_dir: Path):
    """entity stable_id → list of proxy reaction UUIDs (producing/consuming).

    Curated species (often Complexes) that were expanded into virtual variants
    during generation have no UUID of their own in the primary mapping. The
    generator's entity_reaction_proxy_mapping.csv points each such species at the
    UUIDs of the reaction that produces it, so we can read reaction flux as a
    proxy for the species' state. Absent file → no proxies (older catalogs).
    """
    out = defaultdict(list)
    f = pathway_dir / "entity_reaction_proxy_mapping.csv"
    if not f.exists():
        return out
    with open(f) as fh:
        for row in csv.DictReader(fh):
            out[str(row["entity_stable_id"])].append(str(row["proxy_uuid"]))
    return out


def build_adjacency(pathway_dir: Path) -> dict:
    """Forward adjacency: source → [target,...]. Used for reachability checks."""
    adj = defaultdict(list)
    with open(pathway_dir / "logic_network.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            adj[row["source_id"]].append(row["target_id"])
    return adj


def root_occurrences(uuids: list, indeg: Mapping) -> list:
    """The members of `uuids` with no incoming edge: root inputs that are or
    contain the gene (specs/023, the MP-BioPath protocol). Order kept,
    duplicates dropped. Empty when the gene has no root form."""
    seen = set()
    return [u for u in uuids if indeg.get(u, 0) == 0 and not (u in seen or seen.add(u))]


DERIVED_IN_EDGES = {"depletion", "dissociation"}


def recycled_root_occurrences(uuids: list, own_uuids: set, inc: Mapping, fwd: Mapping) -> list:
    """specs/024. For a gene with NO true root: the resolved nodes fed only by
    their own downstream (a catalytic cycle regenerating them) or by derived
    dissociation/depletion edges. The gene's own entity nodes are preferred;
    otherwise the qualifying complexes that are not downstream of another
    qualifying one (or, if they all regenerate each other, the whole cycle's
    pool). `inc[u]` is [(source, edge_type)], `fwd[u]` is [target]."""
    cache: dict = {}

    def reach(u):
        if u not in cache:
            seen, stack = set(), list(fwd.get(u, ()))
            while stack:
                n = stack.pop()
                if n in seen:
                    continue
                seen.add(n)
                stack.extend(fwd.get(n, ()))
            cache[u] = seen
        return cache[u]

    def fed_only_by_itself(u):
        return all(et in DERIVED_IN_EDGES or s in reach(u) for s, et in inc.get(u, ()))

    uniq = list(dict.fromkeys(uuids))
    own = [u for u in uniq if u in own_uuids and fed_only_by_itself(u)]
    if own:
        return own
    cand = [u for u in uniq if fed_only_by_itself(u)]
    entry = [u for u in cand if not any(u in reach(v) for v in cand if v != u)]
    return entry or cand


def entry_occurrences(uuids: list, adj: dict) -> list:
    """The members of `uuids` that no OTHER member reaches (specs/023).

    A member reachable from another is downstream of the gene's entry into the
    network and should be computed, not pinned. If every member is reached
    (they all sit in one cycle), there is no entry point, so all are returned
    and the old protocol applies."""
    members = set(uuids)
    reached = set()
    for start in members:
        seen = {start}
        stack = list(adj.get(start, ()))
        while stack:
            n = stack.pop()
            if n in seen:
                continue
            seen.add(n)
            if n in members and n != start:
                reached.add(n)
            stack.extend(adj.get(n, ()))
    entry = [u for u in uuids if u not in reached]
    return entry or list(uuids)


def duplicate_activator_edges(edges: list) -> set:
    """Indices of activator edges that repeat a (source, target) already seen.

    Keeps the first occurrence, drops the rest, so the entity contributes once
    to the AND product instead of once per role it plays.
    """
    seen = set()
    drop = set()
    for i, e in enumerate(edges):
        if not e.get("is_positive", True):
            continue
        key = (str(e["parent_uuid"]), str(e["child_uuid"]))
        if key in seen:
            drop.add(i)
        else:
            seen.add(key)
    return drop


def node_identity_stids(pathway_dir: Path) -> dict:
    """uuid -> the stable ids that ARE this node (not the ones it contains).

    From nodes.csv this is the node's own ``diagram_entity_id``. It must not be
    derived by inverting ``load_stid_to_uuids``: that index also files every
    node under each of its ``member_leaves``, so an inverse ``{u: s}`` keeps
    whichever stid came last in set-iteration order. The flagged pairs then
    changed with PYTHONHASHSEED (PIP3: 38 or 73 targets with two or more
    flagged inhibitors, on one catalog), and so did the DS_SKIP_SELF_INH arm
    of specs/012."""
    out: dict = defaultdict(set)
    nodes_csv = pathway_dir / "nodes.csv"
    if nodes_csv.exists():
        with open(nodes_csv) as f:
            for row in csv.DictReader(f):
                de = (row.get("diagram_entity_id") or "").strip()
                if de:
                    out[str(row["uuid"])].add(de)
        return out
    with open(pathway_dir / "stid_to_uuid_mapping.csv") as f:
        for row in csv.DictReader(f):
            out[str(row["uuid"])].add(str(row["stable_id"]))
    return out


def self_contained_inhibitor_pairs(pathway_dir: Path) -> set:
    """(source, target) of inhibitor edges whose source contains an activator
    of the SAME target node. Uses the containment table the generator ships, so
    nothing is inferred here. Activators are collected per target uuid, not per
    target stid: pooling by stid would flag an inhibitor on one variant because
    it contains the activator of a sibling variant."""
    cf = pathway_dir / "containment.csv"
    if not cf.exists():
        return set()
    ident = node_identity_stids(pathway_dir)
    contains: dict = {}
    with open(cf) as f:
        for row in csv.DictReader(f):
            contains.setdefault(row["stable_id"], set()).add(row["contains_stable_id"])
    acts: dict = defaultdict(set)
    inh_rows = []
    with open(pathway_dir / "logic_network.csv") as f:
        for row in csv.DictReader(f):
            su, tu = str(row["source_id"]), str(row["target_id"])
            if row.get("pos_neg") == "neg":
                inh_rows.append((su, tu))
            else:
                acts[tu] |= ident.get(su, set())
    pairs = set()
    for su, tu in inh_rows:
        inside = set()
        for sid in ident.get(su, ()):
            inside |= contains.get(sid, set()) - {sid}
        if inside & acts.get(tu, set()):
            pairs.add((su, tu))
    return pairs


def sibling_regulator_pairs(pathway_dir: Path) -> set:
    """(source, target) regulator edges that name a set member which only a
    SIBLING variant of the target reaction uses (specs/025).

    For target T (a variant of reaction R), the members that distinguish its
    siblings are the input leaves some other variant of R has and T lacks.
    The edge is dropped when the regulator's leaves hit one of those and hit
    none of T's own distinguishing leaves."""
    nodes = {}
    with open(pathway_dir / "nodes.csv", newline="") as f:
        for r in csv.DictReader(f):
            leaves = {x for x in (r.get("member_leaves") or "").split("|") if x}
            de = (r.get("diagram_entity_id") or "").strip()
            nodes[r["uuid"]] = (de, leaves | ({de} if de else set()))
    rows = list(csv.DictReader(open(pathway_dir / "logic_network.csv", newline="")))
    in_leaves: dict = defaultdict(set)
    for r in rows:
        if r["pos_neg"] == "pos" and r.get("edge_type") in ("input", "catalyst", "assembly"):
            in_leaves[r["target_id"]] |= nodes.get(r["source_id"], ("", set()))[1]
    variants: dict = defaultdict(list)
    for u, (de, _) in nodes.items():
        if de and u in in_leaves:
            variants[de].append(u)
    drop = set()
    for r in rows:
        if r.get("edge_type") != "regulator":
            continue
        s_, t_ = r["source_id"], r["target_id"]
        de = nodes.get(t_, ("", set()))[0]
        sibs = [v for v in variants.get(de, []) if v != t_]
        if not sibs:
            continue
        mine = in_leaves[t_]
        common = set.intersection(mine, *(in_leaves[v] for v in sibs))
        own_specific = mine - common
        sib_specific = set().union(*(in_leaves[v] for v in sibs)) - mine
        reg = nodes.get(s_, ("", set()))[1]
        if reg & sib_specific and not reg & own_specific:
            drop.add((s_, t_))
    return drop


def load_edge_pairs(pathway_dir: Path, edge_types: set) -> set:
    """(source_id, target_id) pairs whose edge_type is in `edge_types`."""
    pairs = set()
    with open(pathway_dir / "logic_network.csv") as f:
        for row in csv.DictReader(f):
            if row.get("edge_type") in edge_types:
                pairs.add((str(row["source_id"]), str(row["target_id"])))
    return pairs


def reachable_from(adj, sources):
    if not sources:
        return set()
    visited = set(sources)
    frontier = list(sources)
    while frontier:
        nxt = []
        for u in frontier:
            for v in adj.get(u, ()):
                if v not in visited:
                    visited.add(v)
                    nxt.append(v)
        frontier = nxt
    return visited


def categorize_failure(predicted, expected, perturbed_uuids, keyoutput_uuids, reachable):
    """Same buckets the naive baseline uses, so failure modes are comparable."""
    if predicted == expected:
        return "pass"
    if not perturbed_uuids:
        return "gene_not_in_network"
    if not keyoutput_uuids:
        return "keyoutput_not_in_network"
    if not any(ko in reachable for ko in keyoutput_uuids):
        return "no_path"
    if expected == NORMAL and predicted != NORMAL:
        return "false_positive_change"
    return "propagator_missed"


def neo4j_gene_to_stids(genes):
    """Resolve gene names → {gene: [stable_ids]} via the local Reactome Neo4j."""
    from py2neo import Graph
    graph = Graph(NEO4J_URL, auth=(NEO4J_USER, NEO4J_PASSWORD))
    rows = graph.run(
        """
        UNWIND $names AS gene
        MATCH (re:ReferenceEntity)<-[:referenceEntity]-(pe:PhysicalEntity)
        WHERE gene IN re.geneName
        RETURN gene AS gene, COLLECT(DISTINCT pe.stId) AS stids
        """,
        names=list(genes),
    ).data()
    return {r["gene"]: r["stids"] for r in rows}


def neo4j_dbid_to_stid(dbids):
    """Resolve numeric dbIds → {dbId: stId} via the local Reactome Neo4j.

    key_output values in the ground-truth files are numeric dbIds. The stId is
    NOT always ``R-HSA-{dbId}`` — small molecules / species-agnostic entities use
    ``R-ALL-`` (e.g. PIP3 is R-ALL-179838), so a hardcoded R-HSA prefix silently
    fails to find nodes that ARE in the network. Resolve the real stId here.
    """
    from py2neo import Graph
    graph = Graph(NEO4J_URL, auth=(NEO4J_USER, NEO4J_PASSWORD))
    ints = []
    for d in dbids:
        try:
            ints.append(int(d))
        except (TypeError, ValueError):
            continue
    if not ints:
        return {}
    rows = graph.run(
        """
        UNWIND $dbids AS dbid
        MATCH (e:DatabaseObject {dbId: dbid})
        RETURN dbid AS dbid, e.stId AS stid
        """,
        dbids=ints,
    ).data()
    return {str(r["dbid"]): r["stid"] for r in rows if r["stid"]}


def neo4j_set_member_stids(stids):
    """For each EntitySet stId, its member species stIds (all levels, matching
    granularity). Used to resolve a set-typed key-output to the member nodes the
    network actually carries — sets no longer survive as nodes now that the
    generator expands them (see LNG _matching_leaves)."""
    from py2neo import Graph
    graph = Graph(NEO4J_URL, auth=(NEO4J_USER, NEO4J_PASSWORD))
    rows = graph.run(
        """
        UNWIND $stids AS s
        MATCH (e {stId: s}) WHERE e:EntitySet
        MATCH (e)-[:hasMember|hasCandidate*1..]->(m)
        RETURN s AS setid, collect(DISTINCT m.stId) AS members
        """,
        stids=list(stids),
    ).data()
    return {r["setid"]: r["members"] for r in rows if r["members"]}


def neo4j_gtpase_readout_stids(stids):
    """For network stIds, find the pure GTP-bound (active) / GDP-bound (inactive)
    GTPase forms, keyed by the GTPase gene symbol.

    Returns {"GTP": {gene: [stid, ...]}, "GDP": {gene: [stid, ...]}}. A node
    counts only if its displayName is exactly ``GENE:GTP [compartment]`` (or
    ``:GDP``) — i.e. the bare active/inactive form, not a downstream
    GTPase:effector complex — so the readout is the cycle's actual output.
    """
    from py2neo import Graph
    import re
    graph = Graph(NEO4J_URL, auth=(NEO4J_USER, NEO4J_PASSWORD))
    rows = graph.run(
        """
        UNWIND $stids AS s
        MATCH (e:DatabaseObject {stId: s})
        WHERE e.displayName CONTAINS ':GTP' OR e.displayName CONTAINS ':GDP'
        RETURN e.stId AS stid, e.displayName AS name
        """,
        stids=list(stids),
    ).data()
    out = {"GTP": {}, "GDP": {}}
    pat = re.compile(r"^([A-Z0-9]+):(GTP|GDP)(?: \[|$)")
    for r in rows:
        m = pat.match(r["name"] or "")
        if not m:
            continue
        gene, tag = m.group(1), m.group(2)
        out[tag].setdefault(gene, []).append(r["stid"])
    return out


def parse_perturbation_columns(header):
    """Curator columns like 'KRAS_0' / 'KRAS_2' → list of (gene, direction, col)."""
    out = []
    for col in header:
        if col in {"key_output", "control"}:
            continue
        m = re.match(r"^(.+)_(0|2)$", col)
        if m:
            out.append((m.group(1), int(m.group(2)), col))
    return out


def load_curator(pathway_name: str, ground_truth: str = "curator"):
    aliases = ("key_output", "key output", "key_outout", "key outout", "key_ouput")
    if ground_truth == "experimental":
        p = EXPERIMENTAL_DIR / f"{pathway_name}_experimental_results.tsv"
    else:
        p = CURATOR_DIR / f"{pathway_name}_reactome_curator_results.tsv"
    if not p.exists():
        return None
    with open(p) as f:
        # Pandas-free: tolerate variable column counts by reading raw and aligning.
        lines = [ln.rstrip("\n").split("\t") for ln in f if ln.strip()]
    if not lines:
        return None
    header = lines[0]
    key_col = None
    for alias in aliases:
        if alias in header:
            key_col = header.index(alias)
            if alias != "key_output":
                header[key_col] = "key_output"
            break
    if key_col is None:
        return None
    rows = []
    ncol = len(header)
    for line in lines[1:]:
        if len(line) < ncol:
            line = line + [""] * (ncol - len(line))
        rows.append(dict(zip(header, line)))
    return header, rows


def solve_via_ds_api(network_json: dict, observations: dict, network_id=None):
    """Call deltasignal /api/solve.

    When ``network_id`` is given (from a prior /api/parse), send only the
    observations and let the server reuse its cached network — avoids
    re-shipping/re-parsing the whole network on every solve (critical for large
    networks). Falls back to sending the full network otherwise.
    """
    if network_id is not None:
        payload = {"network_id": network_id, "observations": observations}
    else:
        payload = {"network": network_json, "observations": observations}
    body = json.dumps(payload).encode()
    req = Request(f"{DS_BASE}/api/solve",
                  data=body,
                  headers={"Content-Type": "application/json"},
                  method="POST")
    timeout_s = int(os.environ.get("DS_SOLVE_TIMEOUT", "600"))
    with urlopen(req, timeout=timeout_s) as r:
        return json.loads(r.read())


def parse_pathway_via_ds_api(pathway_dir_name: str) -> dict:
    body = json.dumps({"pathway_id": pathway_dir_name}).encode()
    req = Request(f"{DS_BASE}/api/parse",
                  data=body,
                  headers={"Content-Type": "application/json"},
                  method="POST")
    with urlopen(req, timeout=120) as r:
        return json.loads(r.read())


def run_pathway(pathway_id: str, pathway_name: str, gene_to_stids_cache=None,
                ground_truth: str = "curator") -> dict:
    pathway_dir = find_pathway_dir(pathway_id)
    if pathway_dir is None:
        return {"status": "no_network", "name": pathway_name, "id": pathway_id}

    curator = load_curator(pathway_name, ground_truth=ground_truth)
    if curator is None:
        return {"status": f"no_{ground_truth}_file", "name": pathway_name, "id": pathway_id}
    header, rows = curator

    stid_to_uuids = load_stid_to_uuids(pathway_dir)
    entity_reaction_proxies = load_entity_reaction_proxies(pathway_dir)
    perturbations = parse_perturbation_columns(header)

    gene_names = sorted({g for g, _, _ in perturbations})

    if gene_to_stids_cache is not None and pathway_name in gene_to_stids_cache:
        gene_to_stids = gene_to_stids_cache[pathway_name]
    else:
        gene_to_stids = neo4j_gene_to_stids(gene_names)
        if gene_to_stids_cache is not None:
            gene_to_stids_cache[pathway_name] = gene_to_stids

    gene_to_uuids = {}
    for g in gene_names:
        uuids = []
        for sid in gene_to_stids.get(g, []):
            uuids.extend(stid_to_uuids.get(sid, []))
        gene_to_uuids[g] = uuids
    indeg = Counter()
    with open(pathway_dir / "logic_network.csv") as f:
        for row in csv.DictReader(f):
            indeg[row["target_id"]] += 1
    if PIN_SCOPE == "entry":
        pin_adj = build_adjacency(pathway_dir)
        gene_to_uuids = {g: entry_occurrences(us, pin_adj) if us else us
                         for g, us in gene_to_uuids.items()}
    elif PIN_SCOPE in ("root", "root_cycle"):
        roots = {g: root_occurrences(us, indeg) for g, us in gene_to_uuids.items()}
        if PIN_SCOPE == "root_cycle" and any(us and not roots[g] for g, us in gene_to_uuids.items()):
            inc, fwd, ident = defaultdict(list), defaultdict(list), {}
            with open(pathway_dir / "logic_network.csv") as f:
                for row in csv.DictReader(f):
                    inc[row["target_id"]].append((row["source_id"], row.get("edge_type", "")))
                    fwd[row["source_id"]].append(row["target_id"])
            nodes_csv = pathway_dir / "nodes.csv"
            if nodes_csv.exists():
                with open(nodes_csv) as f:
                    ident = {r["uuid"]: r["diagram_entity_id"] for r in csv.DictReader(f)}
            for g, us in gene_to_uuids.items():
                if us and not roots[g]:
                    own_stids = set(gene_to_stids.get(g, []))
                    own = {u for u in us if ident.get(u) in own_stids}
                    roots[g] = recycled_root_occurrences(us, own, inc, fwd)
        gene_to_uuids = roots
    # Record what was actually pinned, so a run's protocol is on the record
    # rather than inferred from its flags (specs/023: 89% of pins had silently
    # become mid-pathway complexes for two months).
    for us in gene_to_uuids.values():
        PIN_TALLY["perturbations"] += bool(us)
        PIN_TALLY["pinned"] += len(set(us))
        PIN_TALLY["pinned_roots"] += sum(1 for u in set(us) if indeg[u] == 0)

    # Resolve key_output dbIds to their REAL stIds (may be R-ALL-, R-NUL-, not
    # just R-HSA-) so species that are genuinely in the network are found.
    dbid_to_stid = neo4j_dbid_to_stid({str(r["key_output"]) for r in rows})

    key_output_uuids = {}
    for r in rows:
        ko = str(r["key_output"])
        if ko and ko not in key_output_uuids:
            # Try the resolved stId first, then the historical R-HSA-{ko} guess.
            candidate_sids = []
            resolved = dbid_to_stid.get(ko)
            if resolved:
                candidate_sids.append(resolved)
            if f"R-HSA-{ko}" not in candidate_sids:
                candidate_sids.append(f"R-HSA-{ko}")
            uuids = []
            for sid in candidate_sids:
                uuids = stid_to_uuids.get(sid, [])
                if uuids:
                    break
            # Fallback: a curated species expanded into virtual variants has no
            # UUID of its own — read the flux of its producing reaction instead.
            if not uuids:
                for sid in candidate_sids:
                    uuids = entity_reaction_proxies.get(sid, [])
                    if uuids:
                        break
            # Set-typed readout: the generator expands EntitySets to their member
            # species, so the set itself has no node. Aggregate the member
            # species that ARE present in the network.
            if not uuids:
                members = neo4j_set_member_stids(candidate_sids)
                seen_m = set()
                for sid in candidate_sids:
                    for m in members.get(sid, []):
                        if m in seen_m:
                            continue
                        seen_m.add(m)
                        uuids = uuids + stid_to_uuids.get(m, [])
            key_output_uuids[ko] = uuids

    # Optional key-output remap for recurated pathways (see DS_KO_REMAP). Build
    # the perturbed-gene → active/inactive-form lookup from the network's own
    # nodes, plus a pooled fallback for upstream-regulator perturbations.
    gtpase_readout = None
    if KO_REMAP and any(ko in GENERIC_GTPASE_READOUTS for ko in key_output_uuids):
        readout_stids = neo4j_gtpase_readout_stids(list(stid_to_uuids.keys()))
        per_gene = {"GTP": {}, "GDP": {}}
        pooled = {"GTP": [], "GDP": []}
        for tag in ("GTP", "GDP"):
            for g, sids in readout_stids[tag].items():
                us = [u for sid in sids for u in stid_to_uuids.get(sid, [])]
                per_gene[tag][g] = us
                pooled[tag].extend(us)
        gtpase_readout = (per_gene, pooled)
        print(f"    [KO_REMAP] {pathway_name}: mapped generic RhoGTPase readouts "
              f"to {len(per_gene['GTP'])} GTP / {len(per_gene['GDP'])} GDP paralog "
              f"forms (pooled fallback {len(pooled['GTP'])}/{len(pooled['GDP'])} uuids)",
              flush=True)

    parsed = parse_pathway_via_ds_api(pathway_dir.name)
    if parsed.get("status") != "success":
        return {"status": "parse_failed", "name": pathway_name, "error": parsed.get("message")}
    # The server resolves `pathway_id` against ITS OWN catalog mount, while
    # every uuid below comes from DS_CATALOG_ROOT. Nothing checks that those two
    # point at the same catalog build, and uuids are regenerated per build — two
    # builds of the same pathway share NONE. When they diverge, every
    # observation silently matches no node, nothing is perturbed, and every
    # readout reports baseline: a full run of 23,908 cases with a single
    # distinct predicted value, scored and reported as if it meant something.
    # Cost two benchmark runs before it was spotted. Fail loudly instead.
    server_uuids = {str(n["uuid"]) for n in parsed["nodes"]}
    local_uuids = {u for uus in stid_to_uuids.values() for u in uus}
    if local_uuids and not (server_uuids & local_uuids):
        raise SystemExit(
            f"\nCATALOG MISMATCH on {pathway_dir.name}:\n"
            f"  the API server's network has {len(server_uuids)} uuids, this "
            f"script resolved {len(local_uuids)} from\n"
            f"  DS_CATALOG_ROOT={CATALOG_ROOT}, and they share NONE.\n"
            f"  The server loads pathway_id from its own mount — set "
            f"PATHWAY_CATALOG to the same\n"
            f"  catalog and recreate the container, or every prediction will "
            f"be baseline.\n")

    edges = parsed["edges"]
    # Diagnostic A/B: drop synthetic boundary edges (assembly/dissociation) of a
    # given type at solve time, to isolate their effect without regenerating.
    skip_pairs = set()
    if SKIP_EDGE_TYPES:
        skip_pairs |= load_edge_pairs(pathway_dir, SKIP_EDGE_TYPES)
    if SKIP_SELF_CONTAINED_INHIBITORS:
        skip_pairs |= self_contained_inhibitor_pairs(pathway_dir)
    if SIBLING_REGULATORS == "member_only":
        skip_pairs |= sibling_regulator_pairs(pathway_dir)
    if skip_pairs:
        edges = [e for e in edges
                 if (str(e["parent_uuid"]), str(e["child_uuid"])) not in skip_pairs]
    if DEDUP_ACTIVATORS:
        drop = duplicate_activator_edges(edges)
        if drop:
            edges = [e for i, e in enumerate(edges) if i not in drop]
    # containment and cofactor_stids travel with a POSTed network, or an arm that
    # drops edges would also silently switch off the self-inhibitor rule and the
    # bundle's cofactor list (the solve would report self_inhibitor_rule inert).
    network_payload = {"nodes": parsed["nodes"], "edges": edges, "pathways": parsed["pathways"],
                       "containment": parsed.get("containment"),
                       "cofactor_stids": parsed.get("cofactor_stids")}
    # Use the server-cached network by id (fast: send only observations per
    # solve). But if SKIP_EDGE_TYPES modified the edges above, the cached
    # network is stale, so send the full modified payload instead.
    modified = (bool(SKIP_EDGE_TYPES) or SKIP_SELF_CONTAINED_INHIBITORS or DEDUP_ACTIVATORS
                or SIBLING_REGULATORS != "off")
    solve_network_id = parsed.get("network_id") if not modified else None

    total = 0
    correct = 0
    valid_total = 0
    valid_correct = 0
    confusion = Counter()
    case_log = []
    failure_categories = Counter()
    adj = build_adjacency(pathway_dir)
    reachable_cache: dict = {}

    # One solve per (gene, direction). Cache the result then read every
    # key_output for that perturbation.
    for gene, direction, col in perturbations:
        uuids = gene_to_uuids.get(gene, [])
        ds_result = None
        if uuids:
            ui_value = PERTURB_UI_DOWN if direction == DOWN else PERTURB_UI_UP
            obs = {u: [ui_value, PIN_CONFIDENCE] for u in uuids}
            try:
                ds_result = solve_via_ds_api(network_payload, obs, network_id=solve_network_id)
            except (HTTPError, URLError) as e:
                return {"status": "solve_failed", "name": pathway_name, "error": str(e)}
            if ds_result.get("status") != "success":
                return {"status": "solve_error", "name": pathway_name,
                        "error": ds_result.get("message")}

        activities = ds_result.get("node_activities", {}) if ds_result else {}

        for r in rows:
            ko = str(r["key_output"])
            ko_uuids = key_output_uuids.get(ko, [])
            # Remap a generic RhoGTPase readout to the perturbed gene's own
            # active/inactive form, falling back to the pooled union of all such
            # forms when the perturbed gene is an upstream regulator (no self form).
            if gtpase_readout is not None and ko in GENERIC_GTPASE_READOUTS:
                tag = GENERIC_GTPASE_READOUTS[ko]
                per_gene, pooled = gtpase_readout
                # Perturbed gene IS a GTPase → read its own active/inactive form
                # (clean, faithful). Upstream regulator (GEF/GAP/GDI) → the pooled
                # union (the original generic-readout semantics). Note: regulator
                # cases are limited not by this mapping but by a genuine coarse→fine
                # divergence — v97 encodes regulatory redundancy (multiple GEFs/GAPs
                # per GTPase, OR logic) that the pooled 2019 curator model lacked, so
                # single-regulator perturbation legitimately doesn't move the readout.
                ko_uuids = per_gene[tag].get(gene) or pooled[tag]
            try:
                expected = int(r[col])
            except (KeyError, ValueError):
                continue
            if expected == -999:
                continue
            # Predict by aggregating the key-output's UUID activities. A species
            # can map to many UUIDs (position-aware variants / multiple producing
            # reactions); how we collapse them matters. DS_KO_AGG selects:
            #   max     — any context active ⇒ present (default; lax for knockouts)
            #   mean    — average across contexts
            #   min     — all contexts must hold (strict; sensitive to knockouts)
            #   extreme — the value deviating most from baseline (handles UP&DOWN)
            if uuids and ko_uuids:
                vals = [activities.get(u, 0.01) for u in ko_uuids]
                if KO_AGG == "mean":
                    agg = sum(vals) / len(vals)
                elif KO_AGG == "min":
                    agg = min(vals)
                elif KO_AGG == "extreme":
                    agg = max(vals, key=lambda v: abs(v - 0.01))
                else:
                    agg = max(vals)
                pred_ui = agg * 100.0
                predicted = classify(pred_ui)
            else:
                pred_ui = 1.0  # default NORMAL
                predicted = NORMAL  # gene or key_output not in network

            total += 1
            is_valid = bool(uuids) and bool(ko_uuids)
            if is_valid:
                valid_total += 1
            if predicted == expected:
                correct += 1
                if is_valid:
                    valid_correct += 1
            confusion[(predicted, expected)] += 1
            cat = "pass"
            if predicted != expected:
                if gene not in reachable_cache:
                    reachable_cache[gene] = reachable_from(adj, set(uuids))
                cat = categorize_failure(predicted, expected, uuids, ko_uuids,
                                         reachable_cache[gene])
                failure_categories[cat] += 1
            case_log.append((gene, direction, ko, predicted, expected, is_valid,
                             len(uuids), len(ko_uuids), cat, pred_ui,
                             "|".join(sorted(set(uuids))), "|".join(sorted(set(ko_uuids)))))

    return {
        "status": "ok",
        "name": pathway_name,
        "id": pathway_id,
        "total": total,
        "correct": correct,
        "accuracy": correct / total if total else 0.0,
        "valid_total": valid_total,
        "valid_correct": valid_correct,
        "valid_accuracy": valid_correct / valid_total if valid_total else 0.0,
        "confusion": dict(confusion),
        "failure_categories": dict(failure_categories),
        "n_perturbations": len(perturbations),
        "n_key_outputs": len(key_output_uuids),
        "case_log": case_log,
    }


def network_edge_count(pathway_id: str) -> int:
    d = find_pathway_dir(pathway_id)
    if d is None:
        return -1
    f = d / "logic_network.csv"
    if not f.exists():
        return -1
    # Subtract 1 for the header row.
    with open(f) as fh:
        return sum(1 for _ in fh) - 1


def summary_metrics(confusion: Mapping[tuple, int]) -> dict:
    """Per-class precision/recall/F1 and the aggregates, from a confusion map.

    `confusion` is keyed (predicted, expected) — the orientation
    `run_pathway` writes. The orientation matters and is easy to get wrong
    silently: F1 is symmetric in precision and recall, so transposing the
    matrix leaves macro-F1 and change-F1 IDENTICAL while turning balanced
    accuracy (macro-recall) into macro-precision. Only the recall-derived
    number moves, which is why it is asserted separately in the tests.

    Extracted from main() so the project's primary decision metric can be
    checked against hand-computed values instead of only ever being observed
    on live data.
    """
    def _prf(cls):
        tp = confusion.get((cls, cls), 0)
        fp = sum(confusion.get((cls, e), 0) for e in (DOWN, NORMAL, UP) if e != cls)
        fn = sum(confusion.get((p, cls), 0) for p in (DOWN, NORMAL, UP) if p != cls)
        prec = tp / (tp + fp) if (tp + fp) else 0.0
        rec = tp / (tp + fn) if (tp + fn) else 0.0
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
        return prec, rec, f1

    prf = {c: _prf(c) for c in (DOWN, NORMAL, UP)}
    f1 = {c: prf[c][2] for c in (DOWN, NORMAL, UP)}
    rec = {c: prf[c][1] for c in (DOWN, NORMAL, UP)}
    return {
        "precision": {c: prf[c][0] for c in (DOWN, NORMAL, UP)},
        "recall": rec,
        "f1": f1,
        "macro_f1": sum(f1.values()) / 3,
        "change_f1": (f1[DOWN] + f1[UP]) / 2,
        "balanced_accuracy": sum(rec.values()) / 3,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pathway-list",
                    default=str(MPBIO_ROOT / "pathway_list.tsv"),
                    help="TSV with pathway_id and pathway_name columns")
    ap.add_argument("--report", default="/tmp/deltasignal_vs_mpbiopath.tsv")
    ap.add_argument("--limit", type=int, default=None,
                    help="Process only the first N pathways (smoke test)")
    ap.add_argument("--max-edges", type=int, default=20000,
                    help="Skip pathways larger than this (default 20k — RAF_MAP's "
                         "1.6M edges chokes the HTTP path)")
    ap.add_argument("--ground-truth", choices=["curator", "experimental"],
                    default="curator",
                    help="Which ground truth to compare against. 'experimental' is the "
                         "smaller higher-quality set (only 10 pathways have it).")
    ap.add_argument("--exclude-experimental", action="store_true",
                    help="Drop the 10 pathways that have experimental ground truth "
                         "(use this for the held-out test after tuning on them).")
    ap.add_argument("--dump-cases", default=None,
                    help="Write a per-case TSV (pathway, gene, direction, key_output, "
                         "predicted, expected, valid, n_gene_uuids, n_ko_uuids, category).")
    args = ap.parse_args()

    excluded = set()
    if args.exclude_experimental:
        # Pathway names that have an experimental_results file
        for p in EXPERIMENTAL_DIR.glob("*_experimental_results.tsv"):
            name = p.name.removesuffix("_experimental_results.tsv")
            excluded.add(name)

    # Read the pathway list (same one MP-BioPath uses).
    with open(args.pathway_list) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx_id = header.index("pathway_id")
        idx_name = header.index("pathway_name")
        pathways = []
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) > idx_id and len(cols) > idx_name:
                pathways.append((cols[idx_id], cols[idx_name]))
    pathways = apply_catalog_ids(pathways, Path(__file__).resolve().parent / "catalog_pathways.tsv")

    # Sort ascending by edge count; skip missing, oversized, and excluded.
    sized = []
    for pid, pname in pathways:
        if pname in excluded:
            print(f"  [skip] {pname:55s} excluded (experimental-set, held-out)", flush=True)
            continue
        n = network_edge_count(pid)
        if n < 0:
            print(f"  [skip] {pname:55s} no network", flush=True)
            continue
        if n > args.max_edges:
            print(f"  [skip] {pname:55s} {n} edges > --max-edges {args.max_edges}", flush=True)
            continue
        sized.append((pid, pname, n))
    sized.sort(key=lambda t: t[2])
    pathways = [(pid, pname) for pid, pname, _ in sized]

    if args.limit:
        pathways = pathways[: args.limit]

    print(f"Running benchmark on {len(pathways)} pathway(s) …", flush=True)
    print(f"Protocol: DS_PIN_SCOPE={PIN_SCOPE} DS_PERTURB_UI_DOWN={PERTURB_UI_DOWN:g} "
          f"DS_PERTURB_UI_UP={PERTURB_UI_UP:g} DS_KO_AGG={KO_AGG} "
          f"DS_SKIP_EDGE_TYPES={','.join(sorted(SKIP_EDGE_TYPES)) or '-'} "
          f"DS_SIBLING_REGULATORS={SIBLING_REGULATORS}", flush=True)
    cache = {}
    results = []
    grand_total = 0
    grand_correct = 0
    grand_valid_total = 0
    grand_valid_correct = 0
    grand_confusion = Counter()

    for pid, pname in pathways:
        t0 = time.time()
        res = run_pathway(pid, pname, gene_to_stids_cache=cache,
                          ground_truth=args.ground_truth)
        elapsed = time.time() - t0
        if res["status"] != "ok":
            print(f"  [skip] {pname:50s} {res['status']}", flush=True)
            results.append(res)
            continue
        print(f"  {pname:50s} {res['correct']:5d}/{res['total']:5d} "
              f"({res['accuracy']:.1%}) | valid {res['valid_correct']:5d}/{res['valid_total']:5d} "
              f"({res['valid_accuracy']:.1%}) | {elapsed:.1f}s",
              flush=True)
        grand_total += res["total"]
        grand_correct += res["correct"]
        grand_valid_total += res["valid_total"]
        grand_valid_correct += res["valid_correct"]
        for k, v in res["confusion"].items():
            grand_confusion[k] += v
        results.append(res)

    # --- Imbalance-aware metrics (PRIMARY) ---------------------------------
    # Accuracy is dominated by the majority NORMAL class, so it rewards
    # predicting "no change". macro-F1 (unweighted mean of per-class F1 over
    # DOWN/NORMAL/UP) weights each class equally, so a lazy all-NORMAL
    # predictor scores F1=0 on both change classes. change-F1 (UP+DOWN only)
    # is the "are we actually detecting regulation" guard. These are the
    # numbers to optimize; accuracy is reported below for continuity.
    metrics = summary_metrics(grand_confusion)
    f1 = metrics["f1"]
    macro_f1 = metrics["macro_f1"]
    change_f1 = metrics["change_f1"]
    bal_acc = metrics["balanced_accuracy"]

    print()
    print("=" * 70)
    print(f"DeltaSignal macro-F1 (DOWN/NORM/UP):  {macro_f1:.4f}   "
          f"[PRIMARY — imbalance-aware]")
    print(f"  change-F1 (UP+DOWN only):           {change_f1:.4f}")
    print(f"  balanced accuracy (macro-recall):   {bal_acc:.4f}")
    print(f"  per-class F1: DOWN={f1[DOWN]:.3f}  NORM={f1[NORMAL]:.3f}  UP={f1[UP]:.3f}")
    print(f"DeltaSignal end-to-end accuracy:  "
          f"{grand_correct}/{grand_total} = "
          f"{grand_correct/grand_total*100 if grand_total else 0:.2f}%")
    print(f"Valid-only (gene+ko in network):  "
          f"{grand_valid_correct}/{grand_valid_total} = "
          f"{grand_valid_correct/grand_valid_total*100 if grand_valid_total else 0:.2f}%")
    print()
    print("Confusion (predicted, expected) → count:")
    for (p, e) in sorted(grand_confusion.keys()):
        print(f"  pred={p} exp={e}: {grand_confusion[(p,e)]}")

    grand_failures: Counter = Counter()
    for r in results:
        for k, v in r.get("failure_categories", {}).items():
            grand_failures[k] += v
    if grand_failures:
        print()
        print("Failure categories (where predicted != expected):")
        for k in ("gene_not_in_network", "keyoutput_not_in_network", "no_path",
                  "false_positive_change", "propagator_missed"):
            if k in grand_failures:
                print(f"  {k}: {grand_failures[k]}")

    # Per-pathway report
    with open(args.report, "w") as f:
        f.write("pathway\tid\tstatus\ttotal\tcorrect\taccuracy\tvalid_total\tvalid_correct\tvalid_accuracy\n")
        for r in results:
            if r["status"] != "ok":
                f.write(f"{r['name']}\t{r.get('id','')}\t{r['status']}\t\t\t\t\t\t\n")
                continue
            f.write(f"{r['name']}\t{r['id']}\t{r['status']}\t"
                    f"{r['total']}\t{r['correct']}\t{r['accuracy']:.6f}\t"
                    f"{r['valid_total']}\t{r['valid_correct']}\t{r['valid_accuracy']:.6f}\n")
    print(f"\nPinned: {PIN_TALLY['pinned']} nodes over {PIN_TALLY['perturbations']} "
          f"perturbations, {PIN_TALLY['pinned_roots']} of them roots (DS_PIN_SCOPE={PIN_SCOPE})")
    print(f"\nPer-pathway report: {args.report}")

    if args.dump_cases:
        with open(args.dump_cases, "w") as f:
            f.write("pathway\tgene\tdirection\tkey_output\tpredicted\texpected\t"
                    "valid\tn_gene_uuids\tn_ko_uuids\tcategory\tpred_ui\t"
                    # The uuids actually pinned and read, so a failure can be
                    # classified on the network without re-deriving resolution
                    # (bench/analysis/case_cyclicity.py expects these).
                    "gene_uuids\toutput_uuids\n")
            for r in results:
                if r["status"] != "ok":
                    continue
                for (gene, direction, ko, pred, exp, valid,
                     ng, nk, cat, pred_ui, gus, kus) in r.get("case_log", []):
                    f.write(f"{r['name']}\t{gene}\t{direction}\t{ko}\t{pred}\t{exp}\t"
                            f"{int(valid)}\t{ng}\t{nk}\t{cat}\t{pred_ui:.6f}\t{gus}\t{kus}\n")
        print(f"Per-case dump: {args.dump_cases}")


if __name__ == "__main__":
    main()
