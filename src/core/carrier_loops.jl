"""
Carrier loops (specs/049): a protein bound and released unchanged, E -> E:X -> E.

The generator's `find_pools` rejects such a cycle as a pool (one core-only
signature, so nothing to interconvert) and lists it in `carrier_loops.csv`:
each loop's free form and the in-edges into that form from inside the loop
(its own release). Iterated, the loop has gain exactly 1 and no anchor, so a
closed one (nothing outside it makes the carrier) drifts with the sweep count:
the 14-3-3 dimer in RAF reads 0.55x at 100 sweeps and 0.0077x at 2,000.

Conservation, with bound forms holding a negligible share of the total (the
convention specs/039 amendment 2 uses for pool intermediates), makes free E
equal total E, which its supply sets, not the loop's turnover. That is the
rule `DS_CYCLE_CARRIERS` applies to an enzyme inside a pool; here it applies to
carrier loops outside pools. Bound forms are not touched.

  DS_CARRIER_LOOP_MODE=off     the loops iterate as before (the default until
                               specs/049 amendment 6).
  DS_CARRIER_LOOP_MODE=closed  free forms of CLOSED loops lose their loop
                               in-edges. In the generator's network those are
                               all their producers, so each reads baseline
                               unless pinned; a silo bridge (off by default) or
                               a boundary edge into it would still feed it.
  DS_CARRIER_LOOP_MODE=closed_supply  as `closed`, and every depletion
                               in-edge into a closed free form is dropped too:
                               with bound forms negligible, binding cannot
                               lower free E any more than release raises it
                               (specs/049 amendment 4, the GRB2 case of P10).
  DS_CARRIER_LOOP_MODE=closed_twin  (default) as `closed_supply`, and each closed free
                               form gains an activating OR in-edge
                               (`carrier_supply`) from every root copy of the
                               same entity (same reactome_id, no in-edge, not
                               itself a listed free form): the gene's own level
                               is its supply (specs/049 amendment 5, AURKA and
                               SMAD7 in P11).
  DS_CARRIER_LOOP_MODE=all     the same as `closed` for every listed loop; a
                               fed free form reads only its producers outside
                               the loop.

`DS_CARRIER_COMPONENT_TWIN=on` (specs/084; default since 2026-10-21; `off` restores the earlier solver) adds, under
`closed_twin`: a closed free form with no starting copy of its own entity (a complex
the pathway only regenerates, such as the ATM-carrying DSB complex in Cell Cycle
Checkpoints) gains AND in-edges (`carrier_component`) from the starting copies of
the maximal entities it contains, so it follows its components. One way only.
Amendment 1: only when the form has no activating in-edge left after its loop
edges are dropped; a form a specialization bridge still feeds is left alone.

`DS_CARRIER_LOOP_OPEN=fed` (specs/077; `off` restores the earlier solver) adds, under any mode but
`off` and `all`: a free form of an OPEN loop loses its loop in-edges when it
keeps an activating producer outside every loop listing it (RBBP8 released by
MMEJ is not new RBBP8; SIRT6 still makes it). A listed form whose every
activating producer is a loop step keeps them: there the "loop" is a scaffold
passing a signal along (the IFN alpha/beta receptor recruiting and releasing
STAT2), and `all` cut it. `fed_single` (specs/077 amendment 1) is `fed` that
keeps a release whose reaction has two or more non-cofactor `input` parents: a
displacement (BAD freeing tBID from BCL-2), where the bound share is what moves.
`fed_single_ub` (specs/088) is `fed_single` that does not count ubiquitin as that
second input: a ligase released after ubiquitinating its target (MDM2:MDM4, CBL) is a
carrier release, not a displacement.
`fed_single` is the default since 2026-10-19 (specs/077, adopted by Adam).

An observation still wins: a pinned free form keeps its pin. A bundle without
`carrier_loops.csv` reports "<mode>: no carrier table" and solves as `off`.
"""

const DS_VALID_CARRIER_LOOP_MODES = Set(["off", "closed", "closed_supply", "closed_twin", "all"])
const DS_VALID_CARRIER_LOOP_OPEN = Set(["off", "fed", "fed_single", "fed_single_ub"])

# specs/088: ubiquitin is a donor co-substrate, never a displacer (the generator's
# `_UBIQUITIN_STIDS`, Release 97)
const UBIQUITIN_STIDS = Set(["R-HSA-68524", "R-HSA-113595", "R-HSA-9660007"])

function carrier_loop_open()::String
    value = get(ENV, "DS_CARRIER_LOOP_OPEN", "fed_single")
    value in DS_VALID_CARRIER_LOOP_OPEN || throw(ArgumentError(
        "DS_CARRIER_LOOP_OPEN must be one of $(join(sort(collect(DS_VALID_CARRIER_LOOP_OPEN)), ", ")); got $(repr(value))"))
    return value
end

const CARRIER_COMPONENT_EDGE = "carrier_component"
const DS_VALID_CARRIER_COMPONENT_TWIN = Set(["off", "on"])

"""specs/084: `DS_CARRIER_COMPONENT_TWIN` (default "on" since 2026-10-21, adopted by Adam; "off" restores the earlier solver)."""
function carrier_component_twin()::Bool
    value = get(ENV, "DS_CARRIER_COMPONENT_TWIN", "on")
    value in DS_VALID_CARRIER_COMPONENT_TWIN || throw(ArgumentError(
        "DS_CARRIER_COMPONENT_TWIN must be one of off, on; got $(repr(value))"))
    return value == "on"
end

function carrier_loop_mode()::String
    value = get(ENV, "DS_CARRIER_LOOP_MODE", "closed_twin")
    value in DS_VALID_CARRIER_LOOP_MODES || throw(ArgumentError(
        "DS_CARRIER_LOOP_MODE must be one of $(join(sort(collect(DS_VALID_CARRIER_LOOP_MODES)), ", ")); got $(repr(value))"))
    return value
end

"""
`(network, rule, forms, dropped, twins, open_dropped)`: `network` with the
listed loop in-edges removed under `mode`, the rule to report, how many
distinct free forms lost an edge, how many edges were removed, how many
root-copy supply edges were added (`closed_twin` only), and how many of the
removed edges were open-loop releases (`open = "fed"`, specs/077). Under
`off`, or with no table, the network is returned as given.
"""
function apply_carrier_loops(network::ReactionNetwork, mode::String; open::String = carrier_loop_open())
    mode == "off" && return network, "off", 0, 0, 0, 0
    network.carrier_loops === nothing && return network, "$mode: no carrier table", 0, 0, 0, 0
    drop = Set{Tuple{String, String}}()
    undrained = Set{String}()   # closed_supply: closed free forms that lose their depletion in-edges
    for (free, src, closed) in network.carrier_loops
        (mode == "all" || closed) && push!(drop, (src, free))
        mode in ("closed_supply", "closed_twin") && closed && push!(undrained, free)
    end
    open_drop = Set{Tuple{String, String}}()
    if open in ("fed", "fed_single", "fed_single_ub") && mode != "all"
        # specs/077: an open free form keeps an activating producer that no loop lists
        listed = Dict{String, Set{String}}()
        for (free, src, _) in network.carrier_loops
            push!(get!(listed, free, Set{String}()), src)
        end
        outside = Set{String}()
        for e in network.edges
            e.is_positive && haskey(listed, e.child_uuid) && !(e.parent_uuid in listed[e.child_uuid]) &&
                push!(outside, e.child_uuid)
        end
        # fed_single: a release driven by a second non-cofactor input is a displacement
        displaced = Set{String}()
        if open in ("fed_single", "fed_single_ub")
            cof = cofactor_uuids(network)
            if open == "fed_single_ub"
                for (u, nd) in network.nodes
                    nd.reactome_id !== nothing && nd.reactome_id in UBIQUITIN_STIDS && push!(cof, u)
                end
            end
            ins = Dict{String, Set{String}}()
            for e in network.edges
                e.is_positive && e.edge_type == "input" && !(e.parent_uuid in cof) &&
                    push!(get!(ins, e.child_uuid, Set{String}()), e.parent_uuid)
            end
            for (k, v) in ins
                length(v) >= 2 && push!(displaced, k)
            end
        end
        for (free, src, closed) in network.carrier_loops
            !closed && free in outside && !(src in displaced) && !((src, free) in drop) &&
                push!(open_drop, (src, free))
        end
        union!(drop, open_drop)
    end
    isempty(drop) && isempty(undrained) && return network, mode, 0, 0, 0, 0
    kept = LogicNetworkEdge[]
    hit = Set{String}()
    for e in network.edges
        # only the activating edge the generator listed; a depletion or an
        # inhibition between the same pair is not supply and stays, except that
        # closed_supply drops every depletion into a closed free form
        if (e.is_positive && (e.parent_uuid, e.child_uuid) in drop) ||
           (!e.is_positive && e.edge_type == "depletion" && e.child_uuid in undrained)
            push!(hit, e.child_uuid)
        else
            push!(kept, e)
        end
    end
    n = length(network.edges) - length(kept)
    nopen = count(e -> e.is_positive && (e.parent_uuid, e.child_uuid) in open_drop, network.edges)
    twins = 0
    if mode == "closed_twin"
        # roots of the network as given; a free form's twin is never another free form
        listed = Set(f for (f, _, _) in network.carrier_loops)
        hasin = Set(e.child_uuid for e in network.edges)
        roots = Dict{String, Vector{String}}()
        # an empty id (an inline network can send "") names no entity, so it
        # must not make every id-less root a twin of every id-less free form
        entity(nd) = (nd === nothing || nd.reactome_id === nothing || isempty(nd.reactome_id)) ?
                     nothing : nd.reactome_id
        for (u, nd) in network.nodes
            (entity(nd) === nothing || u in hasin || u in listed) && continue
            push!(get!(roots, entity(nd), String[]), u)
        end
        component = carrier_component_twin()   # specs/084
        cof = isempty(network.cofactor_stids) ? COFACTOR_STIDS : network.cofactor_stids
        # amendment 1: only a form left with no activating in-edge after the drop
        # (a specialization bridge, specs/067, is still a supply)
        supplied = Set(e.child_uuid for e in kept if e.is_positive && e.edge_type != "depletion")
        for free in sort!(collect(undrained))
            rid = entity(get(network.nodes, free, nothing))
            rid === nothing && continue
            for r in sort(get(roots, rid, String[]))
                push!(kept, LogicNetworkEdge(r, free, false, true, 1.0, "carrier_supply"))
                twins += 1
            end
            # specs/084: no starting copy of its own entity, so the closed form (a
            # complex the pathway only regenerates) follows the starting copies of
            # the maximal entities it contains, ANDed as a complex follows its parts
            (component && isempty(get(roots, rid, String[])) && !(free in supplied)) || continue
            inside = Set(c for c in get(network.containment, rid, Set{String}())
                         if c != rid && !(c in cof) && !isempty(get(roots, c, String[])))
            for c in sort!(collect(inside))
                any(o -> o != c && c in get(network.containment, o, Set{String}()), inside) && continue
                for r in sort(roots[c])
                    push!(kept, LogicNetworkEdge(r, free, true, true, 1.0, CARRIER_COMPONENT_EDGE))
                end
            end
        end
    end
    n == 0 && twins == 0 && length(kept) == length(network.edges) - n && return network, mode, 0, 0, 0, 0
    out = ReactionNetwork(network.nodes, kept, network.set_mappings, network.cofactor_stids,
                          network.containment, network.drug_stids, network.pools,
                          network.pathogen_stids, network.carrier_loops, network.node_leaves)
    return out, mode, length(hit), n, twins, nopen
end

"""
The network the solve actually runs under the active `DS_CARRIER_LOOP_MODE`:
for influence scores, which must differentiate the edges the solve used (the
twin's `carrier_supply` edge, specs/059's `recycled_supply` edges and specs/063's
sequestration depletions), not the
loop edges it dropped. Recycle edges (specs/081) are stripped first, as the
solve does: left in, a starting copy with a recycle in-edge stopped being a
root and gained `recycled_supply` edges the solve never had (review
2026-10-22, L1CAM: 13 edges).
"""
solved_carrier_network(network::ReactionNetwork) =
    first(apply_sequestration(first(apply_recycled_twins(
        first(apply_carrier_loops(first(split_recycle_edges(network)), carrier_loop_mode())),
        recycled_twin_mode())), sequestration_mode()))
