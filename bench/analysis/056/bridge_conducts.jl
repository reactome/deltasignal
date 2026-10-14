# specs/056 check: a context_bridge carries its source's value to the root it feeds.
using DeltaSignal; const DS = DeltaSignal
d = "/ctxp/" * ARGS[1]
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"), nothing)
br = [(e.parent_uuid, e.child_uuid) for e in net.edges if e.edge_type == "context_bridge"]
p = DS.steady_state_params()
for (s, t) in br[1:min(4, end)], v in (0.0, 5.0)
    r = DS.solve_steady_state(net, Dict{String,Tuple{Float64,Float64}}(s => (v, 1.0)), p)
    println("pin source $v -> bridged root reads ", round(100r.node_activities[t], sigdigits=4))
end
