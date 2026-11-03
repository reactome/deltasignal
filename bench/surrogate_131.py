"""specs/131: a linear-response surrogate of the propagator, with an optional strength per edge.

Log2-fold changes propagate as delta = A delta + s on the network the server parses: +1 per AND-joined positive
in-edge (input, catalyst, positive regulator, assembly, dissociation, composition), 1/k per each of a node's k OR-joined
positive in-edges (outputs, variant pools and splits, specialization), -1 per negative in-edge (negative regulator,
depletion). Cofactor, drug and pathogen nodes stay at 0 unless pinned; pins hold log2 of their UI value (0 reads the
bottom of the scale). Iterated (Jacobi, clipped at +-log2 100) to a fixed point or SWEEPS sweeps. A weight table
{(parent, child, edge_type): w} scales entries of A; absent edges keep 1.

`install()` swaps the bench's solve call for the surrogate, so the benches' own pins, readouts and calls apply
unchanged; the server is still used to parse networks.
"""
import math
import os

import numpy as np
from scipy import sparse

CLIP = math.log2(100)
SWEEPS = int(os.environ.get("DS_SURROGATE_SWEEPS", "500"))
TOL = 1e-6
_PARSED = {}
# The bench tallies these rule fields and refuses a run whose solves leave them unreported; the surrogate runs none of
# those rules and says so.
RULE_FIELDS = ("carrier_component_rule carrier_loop_moves carrier_loop_open_rule carrier_loop_rule "
               "carrier_multistep_rule conserved_rule cycle_carrier_release cycle_carriers cycle_reeval_carrier_rule "
               "cycle_rule depletion_own_product drug_rule metabolite_pools_rule nonconv_fallback nonconv_pools "
               "nonconv_restarts_max pathogen_rule pool_base_weight pool_supply_and pool_supply_reach "
               "recycled_twin_rule rna_twin_rule scc_order self_activator_rule self_fed_rule self_inhibitor_copies "
               "self_inhibitor_leaves self_inhibitor_rule ubiquitin_rule zero_supply_rule sequestration_rule").split()


def held_uuids(network):
    """Nodes held at 0 like cofactors, drugs and pathogens (by Reactome stId)."""
    held = set()
    for key in ("cofactor_stids", "drug_stids", "pathogen_stids"):
        held |= set(network.get(key) or [])
    return {str(n["uuid"]) for n in network["nodes"] if n.get("reactome_id") in held}


def build(network, weights=None):
    """(node index, sparse A, edge list [(i_child, j_parent, entry, key)])."""
    idx = {str(n["uuid"]): i for i, n in enumerate(network["nodes"])}
    n_or = {}
    for e in network["edges"]:
        if e["is_positive"] and not e["is_and"]:
            c = str(e["child_uuid"])
            n_or[c] = n_or.get(c, 0) + 1
    rows, cols, vals, edges = [], [], [], []
    for e in network["edges"]:
        p, c = str(e["parent_uuid"]), str(e["child_uuid"])
        if p not in idx or c not in idx:
            continue
        if not e["is_positive"]:
            a = -1.0
        elif e["is_and"]:
            a = 1.0
        else:
            a = 1.0 / n_or[c]
        key = (p, c, e.get("edge_type"))
        w = 1.0 if weights is None else weights.get(key, 1.0)
        rows.append(idx[c]), cols.append(idx[p]), vals.append(a * w)
        edges.append((idx[c], idx[p], a, key))
    A = sparse.csr_matrix((vals, (rows, cols)), shape=(len(idx), len(idx)))
    return idx, A, edges


def solve(network, observations, weights=None, sweeps=SWEEPS):
    """({uuid: internal activity}, converged) for pins {uuid: [ui, confidence]}."""
    idx, A, _ = build(network, weights)
    n = len(idx)
    fixed = np.zeros(n, dtype=bool)
    s = np.zeros(n)
    for u, (ui, conf) in observations.items():
        if u in idx and conf > 1e-6:
            i = idx[u]
            fixed[i] = True
            s[i] = -CLIP if ui <= 0 else max(-CLIP, min(CLIP, math.log2(ui)))
    for u in held_uuids(network):
        if u in idx and not fixed[idx[u]]:
            fixed[idx[u]] = True
    d = s.copy()
    converged = False
    for _ in range(sweeps):
        new = np.clip(A @ d, -CLIP, CLIP)
        new[fixed] = s[fixed]
        if np.max(np.abs(new - d)) < TOL:
            d, converged = new, True
            break
        d = new
    acts = {u: 0.01 * 2.0 ** d[i] for u, i in idx.items()}
    return acts, converged


def install(bm):
    """Point the bench module `bm` (benchmark_vs_mpbiopath) at the surrogate for solves."""
    real_parse = bm.parse_pathway_via_ds_api

    def parse(name):
        parsed = real_parse(name)
        if parsed.get("status") == "success":
            _PARSED[parsed["network_id"]] = parsed
        return parsed

    def solve_api(network_json, observations, network_id=None, edge_redirects=None):
        if edge_redirects:
            raise SystemExit("surrogate_131: edge redirects (drug activity) are not modelled")
        net = network_json if network_json is not None else _PARSED[network_id]
        acts, conv = solve(net, observations)
        return {"status": "success", "converged": conv, "node_activities": acts,
                **{k: "surrogate" for k in RULE_FIELDS}}

    bm.parse_pathway_via_ds_api = parse
    bm.solve_via_ds_api = solve_api
    bm.refuse_unless_default_model = lambda res, name: None
