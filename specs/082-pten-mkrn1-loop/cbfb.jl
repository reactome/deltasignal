using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal
d = "/app/pathway_catalog/R-HSA-8878166"
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
pins, outs = split.(split(strip(read("/scratch/cbfb.tsv", String)), '\t'), '|')
nd = CSV.read(joinpath(d, "nodes.csv"), DataFrame); ent = Dict(r.uuid => string(r.diagram_entity_id) for r in eachrow(nd))
println("pins: ", [ent[p] for p in pins])
b = DS.solve_steady_state(net, Dict(String(p) => (0.0, 1.0) for p in pins))
f(u) = b.node_activities[u] / 0.01
inc = Dict{String, Vector{Any}}(); for e in net.edges; push!(get!(inc, e.child_uuid, Any[]), e); end
u = String(outs[1]); seen = Set{String}()
for step in 1:14
    push!(seen, u)
    es = get(inc, u, [])
    println(rpad(ent[u], 15), " ", round(f(u), digits = 3), String(u) in pins ? " [pinned]" : "", "   in: ",
            join([string(e.is_positive ? "+" : "-", e.edge_type, ":", ent[e.parent_uuid], "=", round(f(e.parent_uuid), digits = 3)) for e in es[1:min(end, 6)]], " "))
    cand = [e for e in es if !(e.parent_uuid in seen)]
    isempty(cand) && break
    # follow what pushes it UP: positive edges with the highest source, negative edges with the lowest
    global u = sort(cand, by = e -> e.is_positive ? -f(e.parent_uuid) : f(e.parent_uuid))[1].parent_uuid
end
