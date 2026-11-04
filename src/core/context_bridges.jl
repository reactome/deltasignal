"""
Context bridges and the solver's structural rules (specs/133).

A context build (`bench/context_build.py`, specs/056 and specs/132) joins a pathway to the pathways that feed it with
`context_bridge` edges into the pathway's root copies. The benchmark reads roots without them, so a gene is pinned at
the same forms with and without context. The solver's structural rules did not: a bridge is an in-edge, so a starting
copy fed by one stopped being a starting copy, and every rule that hangs on starting copies (specs/084 component
twins, specs/059 recycled twins, specs/100 RNA twins, carrier loops) changed. In specs/132 that switched off ATM's
component twins in Cell Cycle Checkpoints.

  DS_CONTEXT_BRIDGE_ORDER=first  (default) every rule sees the bridges, as specs/056 and specs/132 ran.
  DS_CONTEXT_BRIDGE_ORDER=last   the bridges are set aside before every structural rule, as recycle edges are
                                 (specs/081), and added back afterwards, as given (`pos,or`), for the solve.

A network without bridges solves identically under both. The solve reports `context_bridge_order` and
`context_bridges` (bridge edges in the network).
"""

const CONTEXT_BRIDGE_EDGE = "context_bridge"

function context_bridge_order()::String
    value = get(ENV, "DS_CONTEXT_BRIDGE_ORDER", "first")
    value in ("first", "last") ||
        throw(ArgumentError("DS_CONTEXT_BRIDGE_ORDER must be first or last; got $(repr(value))"))
    return value
end

"""`(network without context bridges, the bridges)` under `last`; the network unchanged and no edges under `first`."""
function split_context_bridges(network::ReactionNetwork, order::String)
    held = order == "last" ? [e for e in network.edges if e.edge_type == CONTEXT_BRIDGE_EDGE] : LogicNetworkEdge[]
    isempty(held) && return network, held
    kept = [e for e in network.edges if e.edge_type != CONTEXT_BRIDGE_EDGE]
    out = ReactionNetwork(network.nodes, kept, network.set_mappings, network.cofactor_stids,
                          network.containment, network.drug_stids, network.pools,
                          network.pathogen_stids, network.carrier_loops, network.node_leaves,
                          network.rna_stids)
    return out, held
end

"""`network` with `held` added back."""
function restore_context_bridges(network::ReactionNetwork, held::Vector{LogicNetworkEdge})
    isempty(held) && return network
    return ReactionNetwork(network.nodes, vcat(network.edges, held), network.set_mappings, network.cofactor_stids,
                           network.containment, network.drug_stids, network.pools,
                           network.pathogen_stids, network.carrier_loops, network.node_leaves,
                           network.rna_stids)
end
