using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal
d = "/app/pathway_catalog/R-HSA-453279"
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
obs = Dict(p => (80.0, 1.0) for p in ["9dee8c84-63ff-4903-a7e9-d1335ccae245", "b77963c4-f2be-4fe1-919c-768a6de5c50f"])
nd = CSV.read(joinpath(d, "nodes.csv"), DataFrame); ent = Dict(r.uuid => string(r.diagram_entity_id) for r in eachrow(nd))
res = Dict(o => withenv(() -> DS.solve_steady_state(net, obs), "DS_CARRIER_LOOP_OPEN" => o) for o in ("off", "fed_single"))
base = DS.solve_steady_state(net, Dict{String,Tuple{Float64,Float64}}())
f(o, u) = round(res[o].node_activities[u] / 0.01, digits = 3)
inc = Dict{String,Vector{Any}}(); for e in net.edges; push!(get!(inc, e.child_uuid, Any[]), e); end
function show(u, depth)
    println("  "^depth, ent[u], " ", u[1:8], "  off ", f("off", u), "  fed_single ", f("fed_single", u), "  base ", round(base.node_activities[u]/0.01, digits=3))
    depth >= 3 && return
    for e in get(inc, u, [])
        print("  "^(depth+1), e.is_positive ? "+" : "-", e.edge_type, " ")
        show(e.parent_uuid, depth + 1)
    end
end
for u in ("2c110976", "2ff0e394")
    show(first(k for k in keys(net.nodes) if startswith(k, u)), 0)
end
