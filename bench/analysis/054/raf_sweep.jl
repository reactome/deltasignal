using DeltaSignal; const DS = DeltaSignal
dir = "/catalog/R-HSA-5673001"
net = DS.parse_complete_network(joinpath(dir, "logic_network.csv"), joinpath(dir, "stid_to_uuid_mapping.csv"), nothing)
cases = [split(l, "\t") for l in readlines("/s/raf_cases.tsv")]
c = first(c for c in cases if c[2] == "HRAS" && c[4] == "1268261")
pins = split(c[8], "|"); outs = String.(split(c[9], "|"))
ras = [first(c for c in cases if c[2] == "HRAS" && c[4] == "109783")[9]]
rasouts = String.(split(ras[1], "|"))
p = DS.steady_state_params()
for v in (0.0, 0.5, 0.9, 0.99, 1.0, 1.01, 1.1, 2.0, 80.0)
    r = DS.solve_steady_state(net, Dict{String,Tuple{Float64,Float64}}(String(u) => (v, 1.0) for u in pins), p)
    a = r.node_activities
    m(o) = round(100 * sum(a[u] for u in o) / length(o), sigdigits = 4)
    println(rpad("HRAS $v", 12), " RAS:GTP ", m(rasouts), "   p-MAPK dimers ", m(outs), "   conv ", r.converged, "  iters ", r.iterations)
end
