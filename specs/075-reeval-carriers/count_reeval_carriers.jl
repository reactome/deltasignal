# specs/075: pool-carrier pairs the drives re-evaluate, per pathway, at rest.
using DeltaSignal; const DS = DeltaSignal
root = ENV["BUILD"]; tot = 0
for d in sort(readdir(root))
    isfile(joinpath(root, d, "pools.csv")) || continue
    net = redirect_stdout(devnull) do
        DS.parse_complete_network(joinpath(root, d, "logic_network.csv"), joinpath(root, d, "stid_to_uuid_mapping.csv"), nothing)
    end
    r = redirect_stdout(devnull) do
        DS.solve_steady_state(net, Dict{String,Tuple{Float64,Float64}}(), DS.steady_state_params())
    end
    n = r.diagnostics["cycle_reeval_carriers"]
    global tot += n
    n > 0 && println(d, "\t", n, "\tpools ", r.diagnostics["cycle_pools_solved"])
end
println("TOTAL\t", tot)
