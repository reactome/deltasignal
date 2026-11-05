# specs/136: slots the leaf clause adds (DS_SELF_ACTIVATOR_LEAVES=1), per pathway, counted by the solve itself: one
# solve per pathway with no perturbation, reading `self_activator_leaf_pairs`, so the count is on the network the rule
# runs on (after every structural rule; review 2026-11-05, which found the raw parsed network had been counted).
using DeltaSignal; const DS = DeltaSignal
ENV["DS_SELF_ACTIVATOR_LEAVES"] = "1"
root = ENV["BUILD"]; tot = 0
params = DS.SteadyStateParams(1.0, 0.1, 500, 1e-6, "penalty")
for d in sort(readdir(root))
    isfile(joinpath(root, d, "logic_network.csv")) || continue
    net, res = redirect_stdout(devnull) do
        n = DS.parse_complete_network(joinpath(root, d, "logic_network.csv"), joinpath(root, d, "stid_to_uuid_mapping.csv"), nothing)
        n, DS.solve_steady_state(n, Dict{String, Tuple{Float64, Float64}}(), params)
    end
    n = res.diagnostics["self_activator_leaf_pairs"]
    global tot += n
    n > 0 && println(d, "\t", n)
end
println("LEAF_SLOTS\t", tot)
