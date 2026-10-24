using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal
d = "/app/pathway_catalog/R-HSA-1257604"
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
pins, outs = split.(split(strip(read("/scratch/pten.tsv", String)), '\t'), '|')
nd = CSV.read(joinpath(d, "nodes.csv"), DataFrame); ent = Dict(r.uuid => string(r.diagram_entity_id) for r in eachrow(nd))
b = DS.solve_steady_state(net, Dict(String(p) => (80.0, 1.0) for p in pins))
f(u) = round(b.node_activities[u] / 0.01, digits = 4)
inc = Dict{String, Vector{Any}}(); for e in net.edges; push!(get!(inc, e.child_uuid, Any[]), e); end
pip3 = [u for u in keys(net.nodes) if ent[u] == "R-ALL-179838"]
println("PIP3 copies: ", [(u[1:8], f(u)) for u in pip3])
for u in pip3, e in get(inc, u, [])
    e.is_positive && continue
    println("  PIP3 ", u[1:8], " <- ", e.edge_type, " ", ent[e.parent_uuid], " ", e.parent_uuid[1:8], " = ", f(e.parent_uuid))
    # one level up the depleter
    for e2 in get(inc, e.parent_uuid, [])
        println("        <- ", e2.is_positive ? "+" : "-", e2.edge_type, " ", ent[e2.parent_uuid], " ", e2.parent_uuid[1:8], " = ", f(e2.parent_uuid),
                String(e2.parent_uuid) in pins ? " [pinned]" : "")
    end
end
