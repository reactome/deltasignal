# specs/054 T5: per-entity response to 0.1%, 1% and 10% HRAS pushes.
using DeltaSignal; const DS = DeltaSignal
dir = "/catalog/R-HSA-5673001"
net = DS.parse_complete_network(joinpath(dir, "logic_network.csv"), joinpath(dir, "stid_to_uuid_mapping.csv"), nothing)
cases = [split(l, "\t") for l in readlines("/s/raf_cases.tsv")]
c = first(c for c in cases if c[2] == "HRAS" && c[4] == "1268261")
pins = split(c[8], "|")
p = DS.steady_state_params()
sol(v) = DS.solve_steady_state(net, Dict{String,Tuple{Float64,Float64}}(String(u) => (v, 1.0) for u in pins), p).node_activities
b = sol(1.0)
lf(a, u) = log((a[u] + 1e-12) / (b[u] + 1e-12))
sid(u) = something(net.nodes[u].reactome_id, "?")
res = Dict(v => sol(v) for v in (0.99, 0.999, 0.9, 1.1, 1.01))
# per-entity (stid) max |log fold|
for v in (0.999, 0.99, 0.9, 1.01, 1.1)
    a = res[v]; agg = Dict{String,Float64}()
    for u in keys(b); s = first(split(sid(u), "::")); agg[s] = max(get(agg, s, 0.0), abs(lf(a, u))); end
    moved = count(>(1e-3), values(agg)); big = count(>(log(1.15)), values(agg))
    println("HRAS $v: entities moved >0.1% ", moved, ", >15% ", big)
end
# where does 0.9 differ from 10x the 0.99 response? entities by ratio
a9 = res[0.9]; a99 = res[0.99]
agg9 = Dict{String,Float64}(); agg99 = Dict{String,Float64}()
for u in keys(b); s = first(split(sid(u), "::"))
    agg9[s] = max(get(agg9, s, 0.0), abs(lf(a9, u))); agg99[s] = max(get(agg99, s, 0.0), abs(lf(a99, u))); end
println("\nentity  |lf| at 0.99   |lf| at 0.9   ratio (linear ~10)")
for s in sort(collect(keys(agg9)), by = s -> -agg9[s])[1:60]
    r = agg9[s] / max(agg99[s], 1e-15)
    println(rpad(s, 16), lpad(round(agg99[s]; sigdigits = 3), 12), lpad(round(agg9[s]; sigdigits = 3), 12), lpad(round(r; sigdigits = 3), 12))
end
