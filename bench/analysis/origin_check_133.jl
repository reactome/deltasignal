# specs/133 amendment 1: is dropping cross-pathway rule edges the same as running the structural rules on each
# pathway's part alone? For every pathway directory of a context bundle, compares the edge set the rules give
# (solved_carrier_network's rules: material loops, carrier loops, recycled twins, RNA twins, sequestration) on the
# whole bridge-free network after the drop, with the union of the rules' edge sets on each pathway's part.
# Usage (in the dev container): julia --project=/app origin_check_133.jl BUNDLE_DIR
using DeltaSignal
const D = DeltaSignal

rules(n) = first(D.apply_sequestration(first(D.apply_rna_twins(first(D.apply_recycled_twins(
    first(D.apply_carrier_loops(first(D.split_recycle_edges(first(D.apply_material_loops(n)))), D.carrier_loop_mode())),
    D.recycled_twin_mode())), D.rna_twin_mode())), D.sequestration_mode()))
key(e) = (e.parent_uuid, e.child_uuid, e.edge_type, e.is_positive, e.is_and)

function part(n, origin, q)
    keep = Set(u for (u, p) in origin if p == q)
    loops = n.carrier_loops === nothing ? nothing : [l for l in n.carrier_loops if l[1] in keep]
    pools = n.pools === nothing ? nothing : [p for p in n.pools if all(u -> u in keep, p.forms)]
    D.ReactionNetwork(Dict(u => nd for (u, nd) in n.nodes if u in keep),
                      [e for e in n.edges if e.parent_uuid in keep && e.child_uuid in keep],
                      n.set_mappings, n.cofactor_stids, n.containment, n.drug_stids, pools, n.pathogen_stids, loops,
                      n.node_leaves, n.rna_stids, n.leaf_refs)
end

ENV["DS_CONTEXT_BRIDGE_ORDER"] = "last"
bundle = ARGS[1]
bad = 0
for p in sort(readdir(bundle))
    d = joinpath(bundle, p)
    isfile(joinpath(d, "context_origin.csv")) || continue
    net = D.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
    origin = net.node_origin
    nb = first(D.split_context_bridges(net, "last"))
    length(unique(values(origin))) == 1 && continue
    whole = Set(key.(first(D.drop_cross_origin_edges(rules(nb), origin)).edges))
    parts = Set{Any}()
    for q in unique(values(origin))
        union!(parts, key.(rules(part(nb, origin, q)).edges))
    end
    crossed = last(D.drop_cross_origin_edges(rules(nb), origin))
    same = whole == parts
    global bad += !same
    println(p, "\tpathways=", length(unique(values(origin))), "\tcross_dropped=", crossed, "\tidentical=", same,
            same ? "" : "\tonly_whole=$(length(setdiff(whole, parts)))\tonly_parts=$(length(setdiff(parts, whole)))")
end
println("pathways differing: ", bad)
