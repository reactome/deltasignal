using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal
d = "/app/pathway_catalog/R-HSA-1257604"
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
pins, outs = split.(split(strip(read("/scratch/pten.tsv", String)), '\t'), '|')
nd = CSV.read(joinpath(d, "nodes.csv"), DataFrame); ent = Dict(r.uuid => string(r.diagram_entity_id) for r in eachrow(nd))
for (tag, v) in (("rest", 1.0), ("PTEN 80", 80.0))
    b = DS.solve_steady_state(net, Dict(String(p) => (v, 1.0) for p in pins))
    f(u) = round(b.node_activities[u] / 0.01, digits = 4)
    inc = Dict{String, Vector{Any}}(); for e in net.edges; push!(get!(inc, e.child_uuid, Any[]), e); end
    println("== ", tag)
    frontier = [first(u for u in keys(net.nodes) if startswith(u, "9cc87071"))]; seen = Set{String}()
    for depth in 1:4
        nxt = String[]
        for u in frontier
            u in seen && continue; push!(seen, u)
            for e in get(inc, u, [])
                println("  "^depth, ent[u], "=", f(u), " <- ", e.is_positive ? "+" : "-", e.edge_type, " ", ent[e.parent_uuid], " ",
                        e.parent_uuid[1:8], " = ", f(e.parent_uuid), String(e.parent_uuid) in pins ? " [pinned]" : "")
                push!(nxt, e.parent_uuid)
            end
        end
        frontier = nxt
    end
end
