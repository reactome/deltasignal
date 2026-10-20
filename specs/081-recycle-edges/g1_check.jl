using DeltaSignal, CSV, DataFrames, Statistics
const DS = DeltaSignal
old = "/app/pathway_catalog/R-HSA-453279"; new = get(ENV, "NEW_BUNDLE", "/scratch/rcp_g1/R-HSA-453279")
o2s = Dict(r.uuid => r.stable_id for r in eachrow(CSV.read(joinpath(old, "stid_to_uuid_mapping.csv"), DataFrame)))
s2n = Dict{String, Vector{String}}()
for r in eachrow(CSV.read(joinpath(new, "stid_to_uuid_mapping.csv"), DataFrame)); push!(get!(s2n, r.stable_id, String[]), r.uuid); end
netA = DS.parse_complete_network(joinpath(old, "logic_network.csv"), joinpath(old, "stid_to_uuid_mapping.csv"))
netB = DS.parse_complete_network(joinpath(new, "logic_network.csv"), joinpath(new, "stid_to_uuid_mapping.csv"))
println("pools old ", length(netA.pools), " new ", length(netB.pools))
cls(x) = x < 0.85 ? 0 : (x >= 1.15 ? 2 : 1)
tr(us) = unique(vcat([get(s2n, o2s[u], String[]) for u in us]...))
memo = Dict(); tally = Dict{String,Int}(); bygene = Dict{Tuple{String,String},Vector{Int}}()
for l in eachline(ARGS[1])
    ax, gene, dir, ko, expd, pred0, pins, outs = split(l, '\t')
    key = (gene, dir)
    pa = split(pins, '|'); pb = tr(pa)
    if !haskey(memo, key)
        va = Dict(String(p) => (dir == "0" ? 0.0 : 80.0, 1.0) for p in pa)
        vb = Dict(String(p) => (dir == "0" ? 0.0 : 80.0, 1.0) for p in pb)
        memo[key] = (DS.solve_steady_state(netA, va), DS.solve_steady_state(netB, vb))
    end
    a, b = memo[key]
    oa = split(outs, '|'); ob = tr(oa)
    ca = cls(mean(a.node_activities[u] / 0.01 for u in oa)); cb = cls(mean(b.node_activities[u] / 0.01 for u in ob))
    string(ca) != pred0 && (tally["baseline mismatch"] = get(tally, "baseline mismatch", 0) + 1)
    e = parse(Int, expd)
    k = string(ax, ": ", ca == e ? "right" : "wrong", " -> ", cb == e ? "right" : "wrong")
    tally[k] = get(tally, k, 0) + 1
    if ca != cb && ax == "cur"
        v = get!(bygene, key, [0, 0]); (cb == e) ? (v[1] += 1) : (ca == e ? (v[2] += 1) : nothing)
    end
end
println(sort(collect(tally)))
println("changed perturbations (fixed, broken): ", sort(collect(bygene); by = x -> -(x[2][1] + x[2][2])))
b0 = DS.solve_steady_state(netB, Dict{String,Tuple{Float64,Float64}}())
println("new network at rest: pools solved ", b0.diagnostics["cycle_pools_solved"], ", recycle edges ", b0.diagnostics["recycle_edges"],
        ", nodes off baseline >1%: ", count(u -> abs(b0.node_activities[u] / 0.01 - 1) > 0.01, collect(keys(b0.node_activities))))
