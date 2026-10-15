using DeltaSignal; const DS = DeltaSignal
dir = ENV["PW"]
net = DS.parse_complete_network(joinpath(dir, "logic_network.csv"), joinpath(dir, "stid_to_uuid_mapping.csv"), nothing)
pins = split(ENV["PINS"], "|"); v = parse(Float64, ENV["V"])
r = DS.solve_steady_state(net, Dict{String,Tuple{Float64,Float64}}(String(u) => (v, 1.0) for u in pins), DS.steady_state_params())
open(ENV["OUT"], "w") do f
  for (k, x) in r.node_activities; println(f, k, "\t", 100x); end
end
println("conv ", r.converged)
