"""Bridge a curated entity that positional decomposition split in two.

Positional decomposition can emit several uuids for one curated entity. That is
usually right — different occurrences of a shared species genuinely behave
differently — but it sometimes leaves one uuid that RECEIVES signal and has no
outgoing edge, and another that FEEDS reactions and has no incoming one, with
no route between them. The signal arrives at one copy and the reaction reads
the other, so a perturbation has no path at all.

GPVI is the worked example. Its catalyst `R-HSA-442307` is two disconnected
uuids: one with in-degree 0 feeding reaction `442291`, one with in-degree 1
from `SYK:p-VAV` going nowhere. SYK, PIK3CA and PTPN11 therefore cannot reach
the readout `VAV2_Rho/Rac_effectors:GTP` through the catalyst at all.

WHY THIS IS CAPPED, AND WHY EVERY PREVIOUS ATTEMPT FAILED. The shape is
pervasive — 4,415 such entities across 91 of 92 pathways at Release97 — and
closing it wholesale is not a repair, it is a rewiring: bridging all of them
makes 160,090 nodes newly reachable from signals that currently die. The
marginal reach of ONE bridge is wildly skewed:

    min 3   p25 24   median 88   p75 452   p90 2,779   max 39,374

Eight bridges in `Transcriptional_regulation_by_RUNX1` each connect roughly
39,300 of that network's 39,773 nodes. That is hub-flooding, the mode that cost
macro-F1 0.663 -> 0.479 in the cross-pathway stitch, and it is why
one-bridge-per-silo measured -4: it bridged those too. Full merge measured
-204pp and all-pairs produced 2.3M edges.

So the cap is the whole idea. A bridge that reaches a handful of nodes restores
the local route the curation implies; a bridge that reaches the entire network
asserts that one entity's arrival point feeds everything, which the curation
does not say. 38.9% of bridges reach 50 nodes or fewer; the GPVI catalyst
reaches 36.

OFF by default. `DS_SILO_BRIDGE_MAX_REACH=0` disables it entirely, which is the
shipped behaviour; set it to a node count to enable bridges at or below that
marginal reach.
"""

"""Maximum marginal reach a silo bridge may have, or 0 to add none."""
function silo_bridge_max_reach()::Int
    raw = get(ENV, "DS_SILO_BRIDGE_MAX_REACH", "0")
    value = tryparse(Int, raw)
    if value === nothing || value < 0
        throw(ArgumentError(
            "DS_SILO_BRIDGE_MAX_REACH must be a non-negative integer; got $(repr(raw))"))
    end
    return value
end

"""
Edges bridging sink/source pairs whose marginal reach is within the cap.

Returns an empty vector when the feature is off, when a network has no such
pair, or when every candidate reaches too far. Never adds an edge between two
uuids that already have a route.
"""
function silo_bridge_edges(network::ReactionNetwork,
                           max_reach::Int = silo_bridge_max_reach()
                           )::Vector{LogicNetworkEdge}
    out = LogicNetworkEdge[]
    max_reach == 0 && return out

    adjacency = Dict{String, Vector{String}}()
    indegree = Dict{String, Int}()
    outdegree = Dict{String, Int}()
    for edge in network.edges
        push!(get!(adjacency, edge.parent_uuid, String[]), edge.child_uuid)
        outdegree[edge.parent_uuid] = get(outdegree, edge.parent_uuid, 0) + 1
        indegree[edge.child_uuid] = get(indegree, edge.child_uuid, 0) + 1
    end

    by_stid = Dict{String, Vector{String}}()
    for (uuid, node) in network.nodes
        node.reactome_id === nothing && continue
        push!(get!(by_stid, node.reactome_id, String[]), uuid)
    end

    for (_, uuids) in by_stid
        length(uuids) < 2 && continue
        sinks = [u for u in uuids
                 if get(indegree, u, 0) > 0 && get(outdegree, u, 0) == 0]
        sources = [u for u in uuids
                   if get(indegree, u, 0) == 0 && get(outdegree, u, 0) > 0]
        (isempty(sinks) || isempty(sources)) && continue

        sink = argmax(u -> get(indegree, u, 0), sinks)
        before = _reachable(adjacency, sink, max_reach)
        after = _reachable(adjacency, sink, max_reach; extra = sources)
        gain = length(after) - length(before)
        # A gain at the cap may be a truncated walk, so require strictly under.
        gain > max_reach && continue

        for source in sources
            # `or` so the bridge never imposes AND-completeness on the target,
            # and positive because it asserts identity, not regulation.
            push!(out, LogicNetworkEdge(sink, source, false, true, 1.0, "silo_bridge"))
        end
    end
    return out
end

"""Nodes reachable from `start`, stopping once `cap` new nodes are seen."""
function _reachable(adjacency::Dict{String, Vector{String}}, start::String,
                    cap::Int; extra::Vector{String} = String[])::Set{String}
    seen = Set([start])
    stack = String[start]
    for node in extra
        if !(node in seen)
            push!(seen, node)
            push!(stack, node)
        end
    end
    while !isempty(stack)
        current = pop!(stack)
        for next in get(adjacency, current, String[])
            if !(next in seen)
                push!(seen, next)
                push!(stack, next)
                length(seen) > cap + length(extra) + 1 && return seen
            end
        end
    end
    return seen
end
