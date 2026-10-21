using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal
d = "/app/pathway_catalog/R-HSA-1257604"
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
pins, outs = split.(split(strip(read("/scratch/pten.tsv", String)), '\t'), '|')
nd = CSV.read(joinpath(d, "nodes.csv"), DataFrame); ent = Dict(r.uuid => string(r.diagram_entity_id) for r in eachrow(nd))
a = DS.solve_steady_state(net, Dict(String(p) => (1.0, 1.0) for p in pins))
b = DS.solve_steady_state(net, Dict(String(p) => (1.01, 1.0) for p in pins))
f(r, u) = r.node_activities[u] / 0.01
inc = Dict{String, Vector{Any}}(); for e in net.edges; push!(get!(inc, e.child_uuid, Any[]), e); end
start = first(u for u in keys(net.nodes) if ent[u] == "R-HSA-8847980")
u = start; seen = Set{String}()
for step in 1:25
    push!(seen, u)
    println(rpad(ent[u], 15), " ", u[1:8], "  1.0: ", round(f(a, u), digits = 4), "  1.01: ", round(f(b, u), digits = 4), String(u) in pins ? "  [pinned]" : "")
    es = [e for e in get(inc, u, []) if !(e.parent_uuid in seen)]
    isempty(es) && break
    for e in es
        println("      <- ", e.is_positive ? "+" : "-", e.edge_type, " ", ent[e.parent_uuid], " ", e.parent_uuid[1:8], " 1.01: ", round(f(b, e.parent_uuid), digits = 4))
    end
    # follow the in-edge whose source moved most (positive edges: lowest; negative: highest)
    global u = sort(es, by = e -> e.is_positive ? f(b, e.parent_uuid) : -f(b, e.parent_uuid))[1].parent_uuid
end
