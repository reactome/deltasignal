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

  DS_CARRIER_LOOP_MODE=off     (default) the loops iterate as before.
  DS_CARRIER_LOOP_MODE=closed  free forms of CLOSED loops lose their loop
                               in-edges. In the generator's network those are
                               all their producers, so each reads baseline
                               unless pinned; a silo bridge (off by default) or
                               a boundary edge into it would still feed it.
  DS_CARRIER_LOOP_MODE=all     the same for every listed loop; a fed free form
                               reads only its producers outside the loop.

An observation still wins: a pinned free form keeps its pin.
"""

const DS_VALID_CARRIER_LOOP_MODES = Set(["off", "closed", "all"])

function carrier_loop_mode()::String
    value = get(ENV, "DS_CARRIER_LOOP_MODE", "off")
    value in DS_VALID_CARRIER_LOOP_MODES || throw(ArgumentError(
        "DS_CARRIER_LOOP_MODE must be one of $(join(sort(collect(DS_VALID_CARRIER_LOOP_MODES)), ", ")); got $(repr(value))"))
    return value
end

"""
`(network, rule, forms, dropped)`: `network` with the listed loop in-edges
removed under `mode`, the rule to report, how many distinct free forms lost
an edge and how many edges were removed. Under `off`, or with no table, the
network is returned as given.
"""
function apply_carrier_loops(network::ReactionNetwork, mode::String)
    mode == "off" && return network, "off", 0, 0
    network.carrier_loops === nothing && return network, "$mode: no carrier table", 0, 0
    drop = Set{Tuple{String, String}}()
    for (free, src, closed) in network.carrier_loops
        (mode == "all" || closed) && push!(drop, (src, free))
    end
    isempty(drop) && return network, mode, 0, 0
    kept = LogicNetworkEdge[]
    hit = Set{String}()
    for e in network.edges
        # only the activating edge the generator listed; a depletion or an
        # inhibition between the same pair is not supply and stays
        if e.is_positive && (e.parent_uuid, e.child_uuid) in drop
            push!(hit, e.child_uuid)
        else
            push!(kept, e)
        end
    end
    n = length(network.edges) - length(kept)
    n == 0 && return network, mode, 0, 0
    out = ReactionNetwork(network.nodes, kept, network.set_mappings, network.cofactor_stids,
                          network.containment, network.drug_stids, network.pools,
                          network.pathogen_stids, network.carrier_loops)
    return out, mode, length(hit), n
end
