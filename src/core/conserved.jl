"""
Conserved loop participants (specs/035, design D1 of specs/034).

About 9% of reaction-input edges are supplied ONLY by the reaction's own
downstream: every producer of the input node lies in the same cycle. Nothing
outside the loop anchors such a node, so a loop whose AND inputs are all
loop-carried sits on the gain-1 knife-edge (specs/013, 014), and a push can
drain it to the all-zero state. RAF/MAP kinase is the traced case: F-actin,
CNKSR2 and Ca2+ come back only from "Dissociation of RAS:RAF complex", and the
component reads 0 under every RAS/RAF perturbation.

Most of those inputs (78%) are TRANSFORMED inside the loop (RAS:GTP, the RAF
complexes) and carry signal; they are not touched. This rule applies only to
nodes RELEASED UNCHANGED by the loop: every producer consumes a complex that
contains the node's entity. That is a conserved pool, bound and released, like
a cofactor, whose level is set by total abundance rather than by loop flux.

  DS_CONSERVED_MODE=off    (default) previous behaviour.
  DS_CONSERVED_MODE=inert  conserved nodes are pinned at baseline, as cofactors
                           are. An explicit pinning observation still wins.

Different from DS_SCC_BREAK_CATALYST, which froze every recycling catalyst
(including signal-carrying ones) and lost held-out -65.
"""

const DS_VALID_CONSERVED_MODES = Set(["off", "inert"])

function conserved_mode()::String
    value = get(ENV, "DS_CONSERVED_MODE", "off")
    value in DS_VALID_CONSERVED_MODES || throw(ArgumentError(
        "DS_CONSERVED_MODE must be one of $(join(sort(collect(DS_VALID_CONSERVED_MODES)), ", ")); got $(repr(value))"))
    return value
end

"""uuids of conserved loop participants in `network` (see module docstring)."""
function conserved_uuids(network::ReactionNetwork)::Set{String}
    ids = Set{String}(keys(network.nodes))
    for e in network.edges
        push!(ids, e.parent_uuid); push!(ids, e.child_uuid)
    end
    order = sort!(collect(ids))                     # label-independent indexing
    idx = Dict(u => i for (i, u) in enumerate(order))
    adj = [Set{Int}() for _ in order]
    incoming = [LogicNetworkEdge[] for _ in order]
    for e in network.edges
        push!(adj[idx[e.parent_uuid]], idx[e.child_uuid])
        push!(incoming[idx[e.child_uuid]], e)
    end
    comp, _ = tarjan_scc_components(adj)
    size = Dict{Int, Int}()
    for c in comp
        size[c] = get(size, c, 0) + 1
    end
    stid(u) = (n = get(network.nodes, u, nothing); n === nothing ? nothing : n.reactome_id)
    contains(container, s) = container !== nothing && s !== nothing &&
        (container == s || s in get(network.containment, container, Set{String}()))

    out = Set{String}()
    for e in network.edges
        e.edge_type == "input" || continue
        s, t = idx[e.parent_uuid], idx[e.child_uuid]
        (comp[s] == comp[t] && size[comp[s]] > 1) || continue
        su = order[s]
        su in out && continue
        ss = stid(su)
        ss === nothing && continue
        prods = [p for p in incoming[s] if p.edge_type != "depletion"]
        isempty(prods) && continue
        all(p -> comp[idx[p.parent_uuid]] == comp[s], prods) || continue
        released = all(prods) do p
            if p.edge_type == "dissociation"
                return contains(stid(p.parent_uuid), ss)
            end
            # the producing reaction consumes (or is catalysed by) a complex holding ss
            any(q -> q.edge_type in ("input", "catalyst") && contains(stid(q.parent_uuid), ss),
                incoming[idx[p.parent_uuid]])
        end
        released && push!(out, su)
    end
    return out
end
