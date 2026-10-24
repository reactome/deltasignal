"""
Sequestration (specs/063): a complex that goes nowhere lowers its free partner.

Reactome curates many inhibitors as a plain binding: follistatin binds
activin, USP18 binds IFNAR2, IL13RA2 binds IL13, and the complex is used by
nothing. The biology is that the binding takes the partner out of
circulation. The model propagates fold changes and never lowers a partner by
binding it, so such an inhibitor cannot reach anything. This rule adds a
depletion edge from the dead-end complex to its partner.

A DEAD-END BINDING: a reaction node with at least two non-cofactor entity
inputs and exactly one non-cofactor output C, where every out-edge of C goes
into a dissociation sink. A PARTNER: an input with another use (an out-edge
into anything but that reaction or a dissociation sink).

  DS_SEQUESTRATION=off           (default) nothing changes.
  DS_SEQUESTRATION=deadend_sole  only bindings with at least one input used
                                 nowhere else (the sequesterer) and exactly
                                 one partner: C -| partner.
  DS_SEQUESTRATION=deadend       every dead-end binding: C -| each partner.
  DS_SEQUESTRATION=sole_direct   (amendment 1) the bindings of `deadend_sole`,
                                 but A -| partner from each input A used
                                 nowhere else. C is built from the partner,
                                 so C -| partner fed the partner's own fall
                                 back into it (R4: DSB -414 from two edges);
                                 under multiplicative AND the part of C's
                                 change that is not the partner's is A's.

Needs node kinds (`nodes.csv`); a network without them reports
`<mode>: no node table` and is solved as `off`.
"""

const DS_VALID_SEQUESTRATION_MODES = Set(["off", "deadend_sole", "deadend", "sole_direct"])

function sequestration_mode()::String
    value = get(ENV, "DS_SEQUESTRATION", "off")
    value in DS_VALID_SEQUESTRATION_MODES || throw(ArgumentError(
        "DS_SEQUESTRATION must be one of $(join(sort(collect(DS_VALID_SEQUESTRATION_MODES)), ", ")); got $(repr(value))"))
    return value
end

"""
`(network, rule, edges)`: `network` with a depletion edge from each dead-end
binding complex to its partners under `mode`, the rule to report, and how many
edges were added. Under `off` the network is returned as given.
"""
function apply_sequestration(network::ReactionNetwork, mode::String)
    mode == "off" && return network, "off", 0
    network.node_leaves === nothing && return network, "$mode: no node table", 0
    kind(u) = (k = get(network.node_leaves, u, nothing); k === nothing ? "" : first(k))
    cof(u) = (nd = get(network.nodes, u, nothing);
              nd !== nothing && nd.reactome_id !== nothing && nd.reactome_id in network.cofactor_stids)
    out = Dict{String, Vector{LogicNetworkEdge}}()
    inc = Dict{String, Vector{LogicNetworkEdge}}()
    for e in network.edges
        push!(get!(out, e.parent_uuid, LogicNetworkEdge[]), e)
        push!(get!(inc, e.child_uuid, LogicNetworkEdge[]), e)
    end
    sink(u) = kind(u) == "dissociation_sink"
    added = LogicNetworkEdge[]
    have = Set((e.parent_uuid, e.child_uuid) for e in network.edges if e.edge_type == "depletion")
    for r in sort!([u for u in keys(network.nodes) if kind(u) == "reaction"])
        ins = sort!(unique([e.parent_uuid for e in get(inc, r, LogicNetworkEdge[])
                            if e.edge_type == "input" && !cof(e.parent_uuid)]))
        outs = unique([e.child_uuid for e in get(out, r, LogicNetworkEdge[])
                       if e.edge_type == "output" && !cof(e.child_uuid)])
        (length(ins) < 2 || length(outs) != 1) && continue
        c = outs[1]
        kind(c) in ("reaction", "dissociation_sink", "") && continue
        all(e -> sink(e.child_uuid), get(out, c, LogicNetworkEdge[])) || continue
        used(u) = any(e -> e.child_uuid != r && !sink(e.child_uuid), get(out, u, LogicNetworkEdge[]))
        partners = filter(used, ins)
        if mode in ("deadend_sole", "sole_direct")
            (length(partners) == 1 && length(partners) < length(ins)) || continue
        end
        # the depleting nodes: the complex, or under sole_direct the sequesterers
        from = mode == "sole_direct" ? filter(!used, ins) : [c]
        for s in from, p in partners
            (s, p) in have && continue
            push!(have, (s, p))
            push!(added, LogicNetworkEdge(s, p, true, false, 1.0, "depletion"))
        end
    end
    isempty(added) && return network, mode, 0
    net = ReactionNetwork(network.nodes, vcat(network.edges, added), network.set_mappings,
                          network.cofactor_stids, network.containment, network.drug_stids,
                          network.pools, network.pathogen_stids, network.carrier_loops,
                          network.node_leaves, network.rna_stids)
    return net, mode, length(added)
end
