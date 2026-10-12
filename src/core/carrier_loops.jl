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

An observation still wins: a pinned free form keeps its pin. A bundle without
`carrier_loops.csv` reports "<mode>: no carrier table" and solves as `off`.
"""

const DS_VALID_CARRIER_LOOP_MODES = Set(["off", "closed", "closed_supply", "closed_twin", "all"])

function carrier_loop_mode()::String
    value = get(ENV, "DS_CARRIER_LOOP_MODE", "closed_twin")
    value in DS_VALID_CARRIER_LOOP_MODES || throw(ArgumentError(
        "DS_CARRIER_LOOP_MODE must be one of $(join(sort(collect(DS_VALID_CARRIER_LOOP_MODES)), ", ")); got $(repr(value))"))
    return value
end

"""
`(network, rule, forms, dropped, twins)`: `network` with the listed loop
in-edges removed under `mode`, the rule to report, how many distinct free
forms lost an edge, how many edges were removed, and how many root-copy
supply edges were added (`closed_twin` only). Under `off`, or with no table,
the network is returned as given.
"""
function apply_carrier_loops(network::ReactionNetwork, mode::String)
    mode == "off" && return network, "off", 0, 0, 0
    network.carrier_loops === nothing && return network, "$mode: no carrier table", 0, 0, 0
    drop = Set{Tuple{String, String}}()
    undrained = Set{String}()   # closed_supply: closed free forms that lose their depletion in-edges
    for (free, src, closed) in network.carrier_loops
        (mode == "all" || closed) && push!(drop, (src, free))
        mode in ("closed_supply", "closed_twin") && closed && push!(undrained, free)
    end
    isempty(drop) && isempty(undrained) && return network, mode, 0, 0, 0
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
    twins = 0
    if mode == "closed_twin"
        # roots of the network as given; a free form's twin is never another free form
        listed = Set(f for (f, _, _) in network.carrier_loops)
        hasin = Set(e.child_uuid for e in network.edges)
        roots = Dict{String, Vector{String}}()
        for (u, nd) in network.nodes
            (nd.reactome_id === nothing || u in hasin || u in listed) && continue
            push!(get!(roots, nd.reactome_id, String[]), u)
        end
        for free in sort!(collect(undrained))
            nd = get(network.nodes, free, nothing)
            rid = nd === nothing ? nothing : nd.reactome_id
            rid === nothing && continue
            for r in sort(get(roots, rid, String[]))
                push!(kept, LogicNetworkEdge(r, free, false, true, 1.0, "carrier_supply"))
                twins += 1
            end
        end
    end
    n == 0 && twins == 0 && return network, mode, 0, 0, 0
    out = ReactionNetwork(network.nodes, kept, network.set_mappings, network.cofactor_stids,
                          network.containment, network.drug_stids, network.pools,
                          network.pathogen_stids, network.carrier_loops)
    return out, mode, length(hit), n, twins
end
