"""
RNA twins (specs/100): a root copy of an RNA species reads the copies that
transcription (or any reaction) makes of it.

The generator can split one RNA into a ROOT copy (no in-edge) and a PRODUCED copy
(an `output` in-edge). In PIP3 activates AKT signaling, PTEN mRNA's root copy is
what translation reads, while transcription (stimulated by TP53) makes the other
copy, which feeds only the 12 PTEN mRNA:miRNA complexes that repress translation.
TP53 then reaches AKT only through the repressors, with the sign reversed (all 28
curator cases wrong). One mRNA is one pool, whichever reaction made a node copy.

  DS_RNA_TWIN=off  nothing changes (the solver before 2026-10-24).
  DS_RNA_TWIN=on   (default since 2026-10-24, specs/100 adoption) each root copy of an entity listed in the bundle's `rnas.csv`
                   gains an activating OR in-edge (`rna_supply`) from every
                   produced copy of the same stable id. One-way, as
                   `recycled_supply` (specs/059). The edge exists only in the
                   solver, so the benchmark's root detection is unchanged and a
                   pinned root copy keeps its pin.

A bundle without `rnas.csv` reports "on: no RNA table" and adds nothing.
"""

const RNA_SUPPLY_EDGE = "rna_supply"

function rna_twin_mode()::String
    value = get(ENV, "DS_RNA_TWIN", "on")
    value in ("off", "on") || throw(ArgumentError("DS_RNA_TWIN must be off or on; got $(repr(value))"))
    return value
end

"""`(network, rule, edges added)` under `mode`."""
function apply_rna_twins(network::ReactionNetwork, mode::String)
    mode == "off" && return network, "off", 0
    network.rna_stids === nothing && return network, "on: no RNA table", 0
    has_in = Set{String}(); produced = Set{String}()
    for e in network.edges
        push!(has_in, e.child_uuid)
        e.edge_type == "output" && push!(produced, e.child_uuid)
    end
    roots = Dict{String, Vector{String}}(); made = Dict{String, Vector{String}}()
    for u in sort(collect(keys(network.nodes)))
        s = network.nodes[u].reactome_id
        (s === nothing || !(s in network.rna_stids)) && continue
        u in has_in || push!(get!(roots, s, String[]), u)
        u in produced && push!(get!(made, s, String[]), u)
    end
    added = LogicNetworkEdge[]
    for (s, rs) in sort(collect(roots)), r in rs, p in get(made, s, String[])
        push!(added, LogicNetworkEdge(p, r, false, true, 1.0, RNA_SUPPLY_EDGE))
    end
    isempty(added) && return network, "on", 0
    out = ReactionNetwork(network.nodes, vcat(network.edges, added), network.set_mappings,
                          network.cofactor_stids, network.containment, network.drug_stids,
                          network.pools, network.pathogen_stids, network.carrier_loops,
                          network.node_leaves, network.rna_stids)
    return out, "on", length(added)
end
