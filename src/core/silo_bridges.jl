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

MEASURED, AND IT LOSES. On 71 pathways / 17,466 scored curator cases at
Release97, one shared catalog, conditioned on the experiment being unchanged
(it was, in all 17,966 shared cases — this adds edges but perturbs the same
genes):

    cap 0 -> cap 50   14,960 -> 14,887 correct   net -73   (+73 / -146)
                      accuracy 0.8565 -> 0.8523
                      macro-F1 0.8277 -> 0.8240

The mechanism check passed and the aggregate still lost, which is the useful
part. GPVI recovers exactly as predicted (+12): with the bridge, a SYK knockout
drives `VAV2_Rho/Rac_effectors:GTP` to 0.0005 instead of sitting at 1.0000, and
it does so through the catalyst route the curation implies rather than through
GDP depletion. The fix is right about the case that motivated it.

It is wrong about everything else, and the error has one shape: 126 of the 146
losses are NO_CHANGE turning into a change call. The gains have that shape too
(63 of 73), so the bridge is simply a machine for producing change calls, and
roughly two in three are wrong. The cap reduced the flooding — RUNX1 contributes
55 bridges instead of 106, and its hub bridges are excluded — but it did not
invert that ratio. Worst: RHO_GTPases_activate_CIT -22, ROBO -20,
RHO_GTPases_activate_IQGAPs -14. Best: NODAL +16, GPVI +12.

REPLICATED on an independently regenerated catalog (2026-09-15). The first
measurement used a catalog built 2026-07-16, which predates additive diagram
bridges and 21 other commits to the generator. Rebuilt from current main and
re-run:

    stale catalog   14,960 -> 14,887   net -73   (+73 / -146)   GPVI +12
    fresh catalog   18,083 -> 18,006   net -77   (+63 / -140)   GPVI +12

Same sign, same magnitude, same worst offenders (RHO_GTPases_activate_CIT,
ROBO), and GPVI recovers by exactly +12 in both. The conclusion does not depend
on which catalog it was measured on.

This is the fourth silo fix to fail and the most carefully targeted one. Taken
with the other three it is fair to say the positional silo is NOT an accuracy
lever, even when the bridging is restricted to genuinely local repairs. Note
the contrast with cofactor conduction, which won (+37) by REMOVING 53 spurious
change calls; this loses by adding 126. Over-coupling, not under-connection, is
what this benchmark punishes.

OFF by default. `DS_SILO_BRIDGE_MAX_REACH=0` disables it entirely, which is the
shipped behaviour and what the evidence supports; set it to a node count to
enable bridges at or below that marginal reach. It is retained as a measurement
instrument, not as a recommendation. Do not sweep the cap looking for a
positive value — that is tuning a constant on the evaluation set, and the error
shape says the ratio is structural rather than a threshold artifact.
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
