using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal
d = "/app/pathway_catalog/R-HSA-8848021"
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
nd = CSV.read(joinpath(d, "nodes.csv"), DataFrame); ent = Dict(r.uuid => string(r.diagram_entity_id) for r in eachrow(nd))
u = first(x for x in keys(net.nodes) if startswith(x, "ea6e3b29"))
println("in-edges of p-PTK6 ea6e3b29: ", sort(unique([(ent[e.parent_uuid], e.edge_type) for e in net.edges if e.child_uuid == u])))
rel = [x for x in keys(net.nodes) if get(ent, x, "") == "R-HSA-8848873"]
for r in rel
    println("release copy ", r[1:8], " outputs: ", [(ent[e.child_uuid], e.child_uuid[1:8]) for e in net.edges if e.parent_uuid == r])
end
lp = joinpath(d, "carrier_loops.csv")
println("carrier_loops rows: ", isfile(lp) ? countlines(lp) - 1 : "no file")
