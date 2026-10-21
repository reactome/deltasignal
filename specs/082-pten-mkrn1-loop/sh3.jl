using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal
line = first(l for l in eachline("/scratch/rails.tsv") if occursin("SH3GL1", l))
pw, gene, ko, _, pins, outs = split(line, '\t')
d = "/app/pathway_catalog/R-HSA-177929"
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
P = split(pins, '|'); O = split(outs, '|')
nd = CSV.read(joinpath(d, "nodes.csv"), DataFrame); ent = Dict(r.uuid => string(r.diagram_entity_id) for r in eachrow(nd))
println("pins: ", [ent[p] for p in P])
b = DS.solve_steady_state(net, Dict(String(p) => (1.01, 1.0) for p in P))
f(u) = b.node_activities[u] / 0.01
inc = Dict{String, Vector{Any}}(); for e in net.edges; push!(get!(inc, e.child_uuid, Any[]), e); end
u = String(O[1]); seen = Set{String}()
for step in 1:30
    push!(seen, u)
    es = get(inc, u, [])
    println(rpad(ent[u], 15), " ", u[1:6], " ", round(f(u), digits = 4), String(u) in P ? " [pinned]" : "", "   in: ",
            join([string(e.is_positive ? "+" : "-", e.edge_type, ":", ent[e.parent_uuid], "=", round(f(e.parent_uuid), digits = 3)) for e in es[1:min(end, 5)]], " "))
    cand = [e for e in es if !(e.parent_uuid in seen)]
    isempty(cand) && break
    global u = sort(cand, by = e -> e.is_positive ? f(e.parent_uuid) : -f(e.parent_uuid))[1].parent_uuid
end
