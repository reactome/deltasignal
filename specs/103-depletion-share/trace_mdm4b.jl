using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal
d = "/app/pathway_catalog/R-HSA-3700989"
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
sn = DS.solved_carrier_network(net)
cs = CSV.read("/scratch/cases_tp53.tsv", DataFrame; header = [:gene, :ro, :pins, :outs], delim = '\t', types = String)
names = Dict(r.st => r.name for r in eachrow(CSV.read("/scratch/names.tsv", DataFrame; header = [:st, :name], delim = '\t', types = String)))
sid(u) = something(sn.nodes[u].reactome_id, "")
nm(u) = (s = sid(u); string(s, " ", first(get(names, split(s, "::")[1], ""), 60)))
base = DS.solve_steady_state(net, Dict{String, Tuple{Float64, Float64}}())
pins = String.(split(first(cs[cs.gene .== get(ENV, "GENE", "MDM4"), :pins]), '|'))
ko = DS.solve_steady_state(net, Dict(p => (parse(Float64, get(ENV, "VAL", "0.0")), 1.0) for p in pins))
fold(u) = round(get(ko.node_activities, u, NaN) / max(get(base.node_activities, u, NaN), 1e-15), digits = 3)
inn = Dict{String, Vector{Any}}()
for e in sn.edges; push!(get!(inn, e.child_uuid, Any[]), e); end
for u in sort(collect(keys(sn.nodes)), by = nm)
    n = nm(u)
    (occursin("MDM2", n) || occursin("MDM4", n) || occursin(r"^R-HSA-(3222171|69488) ", n)) || continue
    println(n, "  fold ", fold(u))
    for e in get(inn, u, [])
        println("    <- [", e.edge_type, e.is_positive ? "+" : "-", "] ", nm(e.parent_uuid), " fold ", fold(e.parent_uuid))
    end
end
