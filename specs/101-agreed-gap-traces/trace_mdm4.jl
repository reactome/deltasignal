using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal
d = "/app/pathway_catalog/R-HSA-3700989"
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
sn = DS.solved_carrier_network(net)
cs = CSV.read("/scratch/cases.tsv", DataFrame; header = [:gene, :ro, :pins, :outs], delim = '\t', types = String)
names = Dict(r.st => r.name for r in eachrow(CSV.read("/scratch/names.tsv", DataFrame; header = [:st, :name], delim = '\t', types = String)))
sid(u) = something(sn.nodes[u].reactome_id, "")
nm(u) = (s = sid(u); string(s, " ", first(get(names, s, ""), 45)))
base = DS.solve_steady_state(net, Dict{String, Tuple{Float64, Float64}}())
pins = String.(split(first(cs[cs.gene .== "MDM4", :pins]), '|'))
ko = DS.solve_steady_state(net, Dict(p => (0.0, 1.0) for p in pins))
fold(u) = round(get(ko.node_activities, u, NaN) / max(get(base.node_activities, u, NaN), 1e-15), digits = 3)
inn = Dict{String, Vector{Any}}()
for e in sn.edges; push!(get!(inn, e.child_uuid, Any[]), e); end
function show(u, depth, seen)
    (u in seen || depth > 4) && return; push!(seen, u)
    println("  "^depth, nm(u), "  fold ", fold(u), "  base ", round(base.node_activities[u] / 0.01, digits = 3))
    for e in get(inn, u, [])
        println("  "^(depth + 1), "<- [", e.edge_type, e.is_positive ? "+" : "-", "] ", nm(e.parent_uuid), " fold ", fold(e.parent_uuid))
    end
    for e in get(inn, u, []); show(e.parent_uuid, depth + 1, seen); end
end
for u in keys(sn.nodes)
    startswith(sid(u), "R-HSA-6804885::variant") && occursin("6804932", sid(u)) && show(u, 0, Set{String}())
end
