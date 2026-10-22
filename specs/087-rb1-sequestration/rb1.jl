using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal
d = "/app/pathway_catalog/R-HSA-453279"
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
nd = CSV.read(joinpath(d, "nodes.csv"), DataFrame); ent = Dict(r.uuid => string(r.diagram_entity_id) for r in eachrow(nd))
pin = first(u for u in keys(net.nodes) if startswith(u, "09d62d8a"))
r = DS.solve_steady_state(net, Dict(pin => (0.0, 1.0)))
f(u) = round(r.node_activities[u] / 0.01, digits = 4)
for (tag, e) in (("RB1 cytosol", "R-HSA-9660621"), ("RB1 nucleus", "R-HSA-68642"), ("RB1:E2F1-3:DP", "R-HSA-188343"),
                 ("E2F1-3:DP", "R-HSA-1227905"), ("CCNA1", "R-HSA-68891"))
    us = [u for u in keys(net.nodes) if ent[u] == e]
    isempty(us) && (println(tag, " no node"); continue)
    vals = sort(unique([f(u) for u in us]))
    println(rpad(tag, 16), e, " copies ", length(us), " values ", vals[1:min(end, 6)])
end
# what consumes the RB1:E2F products of 9018017?
outs = Set{String}()
for e in net.edges
    startswith(ent[e.parent_uuid], "R-HSA-9018017") && push!(outs, e.child_uuid)
end
println("products of RB1 binding E2F: ", sort(unique([ent[u] for u in outs])), " values ", sort(unique([f(u) for u in outs]))[1:min(end,4)])
users = Dict{String, Set{String}}()
for e in net.edges
    e.parent_uuid in outs || continue
    push!(get!(users, ent[e.parent_uuid], Set{String}()), string(ent[e.child_uuid], ":", e.edge_type, e.is_positive ? "+" : "-"))
end
for (k, v) in users; println("  ", k, " -> ", sort(collect(v))[1:min(end, 8)]); end
println("pools touching RB1: ", [p.id for p in net.pools if any(u -> ent[u] in ("R-HSA-68642", "R-HSA-188343"), vcat(p.forms, p.intermediates))])
