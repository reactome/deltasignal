using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal
d = "/app/pathway_catalog/R-HSA-1257604"
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
nd = CSV.read(joinpath(d, "nodes.csv"), DataFrame); ent = Dict(r.uuid => string(r.diagram_entity_id) for r in eachrow(nd))
pins, outs = split.(split(strip(read("/scratch/pten.tsv", String)), '\t'), '|')
byp(p) = first(u for u in keys(net.nodes) if startswith(u, p))
watch = [("translation", byp("9cc87071")), ("miR complex 2318750", byp("a04983dd")), ("PTEN copy", byp("57a5dff3")),
         ("p-MKRN1", byp("45520a0c")), ("PIP3", byp("d7b2a6f6"))]
for w in ("off", nothing), it in (500, 5000)
    b = withenv("DS_SELF_INHIBITOR_WEIGHT" => w) do
        DS.solve_steady_state(net, Dict(String(p) => (80.0, 1.0) for p in pins), DS.SteadyStateParams(1.0, 0.1, it, 1e-6, "penalty"))
    end
    f(u) = round(b.node_activities[u] / 0.01, digits = 3)
    println("w=", something(w, "0.1"), " iters=", it, " conv=", b.converged, ": ", [(k, f(u)) for (k, u) in watch], " readout ", [f(u) for u in outs])
end
# the complex's producer, under the default
b = DS.solve_steady_state(net, Dict(String(p) => (80.0, 1.0) for p in pins))
f(u) = round(b.node_activities[u] / 0.01, digits = 3)
inc = Dict{String, Vector{Any}}(); for e in net.edges; push!(get!(inc, e.child_uuid, Any[]), e); end
c = byp("035863c9")   # a producer copy of 2318750
println("producer ", ent[c], " = ", f(c))
for e in get(inc, c, []); println("   <- ", e.is_positive ? "+" : "-", e.edge_type, " ", ent[e.parent_uuid], " = ", f(e.parent_uuid), String(e.parent_uuid) in pins ? " [pinned]" : ""); end
