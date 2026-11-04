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

Amendment 1: setting bridges aside is not enough. Rules that join copies of one entity by identity (specs/059 recycled
twins, specs/049 carrier twins, specs/084 component twins, specs/100 RNA twins) also joined a pathway's copy to another
pathway's copy of the same entity: in TP53 transcription, PIP3's imported p-AKT1 gained a twin edge from TP53's own
p-AKT1, a loop back into the context through a shared molecule. Under `last`, a context bundle's `context_origin.csv`
(the pathway each node was built from) is read and every edge a rule drew between two pathways is dropped, so each
pathway's rules see only its own copies. Reported as `context_cross_edges`.
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

"""
Under `last` (`origin` from `context_origin.csv`, else `nothing`): `network` without its edges between nodes of two
pathways, and how many were dropped (specs/133 amendment 1). A context bundle's own edges cross pathways only at
context bridges, which are set aside when this runs, so every edge dropped is one a structural rule drew.
"""
function drop_cross_origin_edges(network::ReactionNetwork, origin::Union{Nothing, Dict{String, String}})
    origin === nothing && return network, 0
    function crosses(e)
        a, b = get(origin, e.parent_uuid, nothing), get(origin, e.child_uuid, nothing)
        return a !== nothing && b !== nothing && a != b
    end
    kept = [e for e in network.edges if !crosses(e)]
    n = length(network.edges) - length(kept)
    n == 0 && return network, 0
    return ReactionNetwork(network.nodes, kept, network.set_mappings, network.cofactor_stids,
                           network.containment, network.drug_stids, network.pools,
                           network.pathogen_stids, network.carrier_loops, network.node_leaves,
                           network.rna_stids), n
end
