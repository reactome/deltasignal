# specs/054 T5: small HRAS pushes at iteration caps 500, 2000, 8000.
using DeltaSignal; const DS = DeltaSignal
dir = "/catalog/R-HSA-5673001"
net = DS.parse_complete_network(joinpath(dir, "logic_network.csv"), joinpath(dir, "stid_to_uuid_mapping.csv"), nothing)
cases = [split(l, "\t") for l in readlines("/s/raf_cases.tsv")]
c = first(c for c in cases if c[2] == "HRAS" && c[4] == "1268261")
pins = split(c[8], "|"); outs = String.(split(c[9], "|"))
ras = [first(c for c in cases if c[2] == "HRAS" && c[4] == "109783")[9]]
rasouts = String.(split(ras[1], "|"))
p = DS.steady_state_params()
for mi in (500, 2000, 8000), v in (0.99, 0.999, 1.001, 1.01)
    p = DS.steady_state_params(max_iters = mi)
    r = DS.solve_steady_state(net, Dict{String,Tuple{Float64,Float64}}(String(u) => (v, 1.0) for u in pins), p)
    a = r.node_activities
    m(o) = round(100 * sum(a[u] for u in o) / length(o), sigdigits = 4)
    println("cap $mi ", rpad("HRAS $v", 12), " RAS:GTP ", m(rasouts), "   p-MAPK dimers ", m(outs), "   conv ", r.converged, "  iters ", r.iterations)
end
