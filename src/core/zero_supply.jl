"""
Zero supply (specs/052): a protein that is not there is in no complex.

Under a knockout the root forms of a protein are pinned to 0, but a positive
cycle can go on regenerating complexes built from it (Signaling by EGFR: the
curated CBL endocytosis cycle rails at 100x with no EGFR in the cell). That
breaks conservation, not a parameter.

For one solve, a leaf L is a CANDIDATE if a node pinned to 0 is built from it.
A node CARRIES L if L is in its `member_leaves` (any kind), or if it is a
reaction node with an `input` edge from an L-carrier. An edge carries L if both
ends carry it and its type moves material (input, output, assembly,
dissociation, set_member, variant_pool, variant_split); catalysts, regulators
and depletion never build their target. An L-SOURCE is a non-reaction carrier
with no carrying in-edge. SUPPLIED = reachable along carrying edges from a
source that is not pinned to 0 (a node pinned above 0 is a supplied source; one
pinned to 0 neither is supplied nor passes supply on).

A node REQUIRES L if it carries it and is an entity built from it
(`simple_entity`, `simple_complex`, `set_variant`, `dissociation_sink`); pools
are OR nodes and reactions are not entities. Every unpinned node that requires
a candidate leaf it is not supplied with is held at 0.

  DS_ZERO_SUPPLY_MODE=off   (default) nothing is held.
  DS_ZERO_SUPPLY_MODE=hold  the rule above.

An observation that pins wins. A bundle without `nodes.csv` reports
"hold: no node table" and solves as `off`.
"""

const DS_VALID_ZERO_SUPPLY_MODES = Set(["off", "hold"])

const ZERO_SUPPLY_CARRYING = Set(["input", "output", "assembly", "dissociation",
                                  "set_member", "variant_pool", "variant_split"])
const ZERO_SUPPLY_REQUIRING = Set(["simple_entity", "simple_complex", "set_variant", "dissociation_sink"])

function zero_supply_mode()::String
    value = get(ENV, "DS_ZERO_SUPPLY_MODE", "off")
    value in DS_VALID_ZERO_SUPPLY_MODES || throw(ArgumentError(
        "DS_ZERO_SUPPLY_MODE must be one of $(join(sort(collect(DS_VALID_ZERO_SUPPLY_MODES)), ", ")); got $(repr(value))"))
    return value
end

"""
`(held, rule, leaves)`: the uuids to hold at 0 under `mode`, the rule to report,
and how many candidate leaves had no supply. `pinned` maps a node to its pinned
UI value (observations that clear the confidence gate).
"""
function zero_supply_uuids(network::ReactionNetwork, pinned::Dict{String, Float64}, mode::String)
    held = Set{String}()
    mode == "off" && return held, "off", 0
    network.node_leaves === nothing && return held, "hold: no node table", 0
    nl = network.node_leaves
    candidates = Set{String}()
    for (u, v) in pinned
        v == 0.0 && haskey(nl, u) && union!(candidates, nl[u][2])
    end
    isempty(candidates) && return held, "hold", 0
    kind(u) = haskey(nl, u) ? nl[u][1] : ""
    isreaction(u) = kind(u) == "reaction"
    leaves(u) = haskey(nl, u) ? nl[u][2] : String[]
    out_edges = Dict{String, Vector{LogicNetworkEdge}}()
    for e in network.edges
        e.edge_type in ZERO_SUPPLY_CARRYING || continue
        push!(get!(out_edges, e.parent_uuid, LogicNetworkEdge[]), e)
    end
    empty_leaves = 0
    for L in sort(collect(candidates))
        # entity carriers, then reactions fed by one through an input edge
        carrier = Set{String}(u for u in keys(nl) if !isreaction(u) && L in leaves(u))
        for u in collect(carrier), e in get(out_edges, u, LogicNetworkEdge[])
            e.edge_type == "input" && isreaction(e.child_uuid) && push!(carrier, e.child_uuid)
        end
        has_in = Set{String}()
        for u in carrier, e in get(out_edges, u, LogicNetworkEdge[])
            e.child_uuid in carrier && push!(has_in, e.child_uuid)
        end
        zero(u) = get(pinned, u, -1.0) == 0.0
        frontier = String[u for u in carrier
                          if !zero(u) && ((!isreaction(u) && !(u in has_in)) || get(pinned, u, 0.0) > 0.0)]
        supplied = Set{String}(frontier)
        while !isempty(frontier)
            u = pop!(frontier)
            for e in get(out_edges, u, LogicNetworkEdge[])
                v = e.child_uuid
                (v in carrier && !(v in supplied) && !zero(v)) || continue
                push!(supplied, v)
                push!(frontier, v)
            end
        end
        lost = String[u for u in carrier
                      if !(u in supplied) && !haskey(pinned, u) && kind(u) in ZERO_SUPPLY_REQUIRING]
        isempty(lost) && continue
        empty_leaves += 1
        union!(held, lost)
    end
    return held, "hold", empty_leaves
end
