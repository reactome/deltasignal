using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal
d = "/app/pathway_catalog/R-HSA-195721"
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
nd = CSV.read(joinpath(d, "nodes.csv"), DataFrame); ent = Dict(r.uuid => string(r.diagram_entity_id) for r in eachrow(nd))
pins = String.(split(ARGS[1], '|')); outs = String.(split(ARGS[2], '|'))
pl = only(p for p in net.pools if p.id == "release1")
println("bases ", length(pl.bases), " forms ", length(pl.forms), " pinned that are pool forms: ", count(in(pl.forms), pins),
        " pinned bases: ", count(in(pl.bases), pins))
for m in ("full", "split"), v in (0.0, nothing)
    obs = v === nothing ? Dict{String,Tuple{Float64,Float64}}() : Dict(p => (v, 1.0) for p in pins)
    r = withenv("DS_POOL_BASE_WEIGHT" => m) do; DS.solve_steady_state(net, obs); end
    f(u) = round(r.node_activities[u] / 0.01, digits = 4)
    nb = [u for u in pl.forms if !(u in pl.bases)]
    println("   split=", get(r.diagnostics, "cycle_bases_split", "?"), " rule=", get(r.diagnostics, "pool_base_weight", "?"), " solved=", get(r.diagnostics, "cycle_pools_solved", "?"))
    println("== $m ", v === nothing ? "rest" : "APC KO", ": readout ", [f(u) for u in outs],
            " non-base forms ", [(ent[u], f(u)) for u in nb],
            " bases unpinned ", sort(unique([f(u) for u in pl.bases if !(u in pins)]))[1:min(end,5)],
            " bases pinned ", count(in(pins), pl.bases))
end
