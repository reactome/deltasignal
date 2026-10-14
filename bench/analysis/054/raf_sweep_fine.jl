# specs/054 T5: HRAS swept finely around baseline (run with and without DS_DEDUP_ACTIVATORS=1).
using DeltaSignal; const DS = DeltaSignal
dir = "/catalog/R-HSA-5673001"
net = DS.parse_complete_network(joinpath(dir, "logic_network.csv"), joinpath(dir, "stid_to_uuid_mapping.csv"), nothing)
cases = [split(l, "\t") for l in readlines("/s/raf_cases.tsv")]
c = first(c for c in cases if c[2] == "HRAS" && c[4] == "1268261")
pins = split(c[8], "|"); outs = String.(split(c[9], "|"))
ras = [first(c for c in cases if c[2] == "HRAS" && c[4] == "109783")[9]]
rasouts = String.(split(ras[1], "|"))
p = DS.steady_state_params()
for v in (0.93, 0.95, 0.96, 0.97, 0.98, 1.02, 1.03, 1.04, 1.05, 1.07)
    r = DS.solve_steady_state(net, Dict{String,Tuple{Float64,Float64}}(String(u) => (v, 1.0) for u in pins), p)
    a = r.node_activities
    m(o) = round(100 * sum(a[u] for u in o) / length(o), sigdigits = 4)
    println(rpad("HRAS $v", 12), " RAS:GTP ", m(rasouts), "   p-MAPK dimers ", m(outs), "   conv ", r.converged, "  iters ", r.iterations)
end
