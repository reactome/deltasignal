using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal
d = "/app/pathway_catalog/R-HSA-1257604"
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
pins, outs = split.(split(strip(read("/scratch/pten.tsv", String)), '\t'), '|')
nd = CSV.read(joinpath(d, "nodes.csv"), DataFrame); ent = Dict(r.uuid => string(r.diagram_entity_id) for r in eachrow(nd))
nm = Dict(r.uuid => string(r.node_kind) for r in eachrow(nd))
println("pins: ", [(ent[p], p[1:8]) for p in pins])
b = DS.solve_steady_state(net, Dict(String(p) => (80.0, 1.0) for p in pins))
println("converged ", b.converged)
f(u) = round(b.node_activities[u] / 0.01, digits = 3)
println("readout copies: ", [(ent[u], f(u)) for u in outs])
inc = Dict{String, Vector{Any}}(); for e in net.edges; push!(get!(inc, e.child_uuid, Any[]), e); end
u = String(outs[1]); seen = Set{String}()
for step in 1:30
    push!(seen, u)
    println(rpad(ent[u], 16), " ", u[1:8], "  ", f(u), String(u) in pins ? "  [pinned]" : "")
    es = [e for e in get(inc, u, []) if !(e.parent_uuid in seen)]
    isempty(es) && break
    for e in es
        println("      <- ", e.is_positive ? "+" : "-", e.edge_type, " ", ent[e.parent_uuid], " ", e.parent_uuid[1:8], " ", f(e.parent_uuid))
    end
    # follow the input that pushes the target up most: highest activator, lowest inhibitor/depleter
    global u = sort(es, by = e -> e.is_positive ? -f(e.parent_uuid) : f(e.parent_uuid))[1].parent_uuid
end
