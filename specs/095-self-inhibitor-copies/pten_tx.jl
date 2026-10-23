using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal
d = "/app/pathway_catalog/R-HSA-1257604"
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
nd = CSV.read(joinpath(d, "nodes.csv"), DataFrame); ent = Dict(r.uuid => string(r.diagram_entity_id) for r in eachrow(nd))
pins, outs = split.(split(strip(read("/scratch/pten.tsv", String)), '\t'), '|')
inc = Dict{String, Vector{Any}}(); for e in net.edges; push!(get!(inc, e.child_uuid, Any[]), e); end
tx = first(u for u in keys(net.nodes) if startswith(u, "551ba5e1"))
for w in ("off", nothing)
    b = withenv("DS_SELF_INHIBITOR_WEIGHT" => w) do; DS.solve_steady_state(net, Dict(String(p) => (80.0, 1.0) for p in pins)); end
    f(u) = round(b.node_activities[u] / 0.01, digits = 3)
    println("== w=", something(w, "0.1"), " ", ent[tx], " = ", f(tx))
    for e in get(inc, tx, []); println("     <- ", e.is_positive ? "+" : "-", e.edge_type, " ", ent[e.parent_uuid], " ", e.parent_uuid[1:8], " = ", f(e.parent_uuid), String(e.parent_uuid) in pins ? " [pinned]" : ""); end
end
