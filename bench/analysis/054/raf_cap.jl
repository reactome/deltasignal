using DeltaSignal; const DS = DeltaSignal
dir = "/catalog/R-HSA-5673001"
net = DS.parse_complete_network(joinpath(dir, "logic_network.csv"), joinpath(dir, "stid_to_uuid_mapping.csv"), nothing)
cases = [split(l, "\t") for l in readlines("/s/raf_cases.tsv")]
c = first(c for c in cases if c[2] == "HRAS" && c[4] == "1268261")
pins = split(c[8], "|"); outs = String.(split(c[9], "|"))
for it in (100, 250, 500, 501, 1000, 2000, 5000)
    p = DS.steady_state_params(max_iters = it)
    for v in (0.0, 80.0)
        r = DS.solve_steady_state(net, Dict{String,Tuple{Float64,Float64}}(String(u) => (v, 1.0) for u in pins), p)
        a = r.node_activities
        println(rpad("cap $it HRAS $v", 22), " p-MAPK ", round(100 * sum(a[u] for u in outs) / length(outs), sigdigits = 4),
                "  conv ", r.converged, "  iters ", r.iterations, "  residual ", get(r.diagnostics, "final_residual", get(r.diagnostics, "residual", "?")))
    end
end
