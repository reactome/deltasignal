# specs/070: one case solved on one bundle, every node's activity written out.
# env PW = bundle dir, PINS = uuids joined by | (empty: no pins), V = pinned
# value (0-100 scale), OUT = file of uuid, activity (0-100 scale).
using DeltaSignal; const DS = DeltaSignal
dir = ENV["PW"]
net = DS.parse_complete_network(joinpath(dir, "logic_network.csv"), joinpath(dir, "stid_to_uuid_mapping.csv"), nothing)
pins = filter(!isempty, split(ENV["PINS"], "|")); v = parse(Float64, ENV["V"])
r = DS.solve_steady_state(net, Dict{String,Tuple{Float64,Float64}}(String(u) => (v, 1.0) for u in pins), DS.steady_state_params())
open(ENV["OUT"], "w") do f
  for (k, x) in r.node_activities; println(f, k, "\t", 100x); end
end
println("conv ", r.converged, " ", join(["$k=$(r.diagnostics[k])" for k in sort(collect(keys(r.diagnostics))) if occursin("nonconv", k) || occursin("cycle_pools", k)], " "))
