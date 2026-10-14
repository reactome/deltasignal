# specs/056 trace: one case in a context bundle. Prints each bridge INTO the
# pathway whose source moved, with the source's pathway and entity, and the
# readouts. ARGS: bundle dir, case row file (key_output, expected, predicted,
# pred_ui, gene_uuids, output_uuids), owner table (uuid, pathway, stid, kind).
using DeltaSignal; const DS = DeltaSignal
d = ARGS[1]
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"), nothing)
own = Dict(l[1] => (l[2], l[3], l[4]) for l in (split(x, "\t") for x in readlines(ARGS[3])))
P = basename(d)
c = split(readlines(ARGS[2])[1], "\t")
pins, outs = String.(split(c[5], "|")), String.(split(c[6], "|"))
p = DS.steady_state_params()
r0 = DS.solve_steady_state(net, Dict{String,Tuple{Float64,Float64}}(), p)
r = DS.solve_steady_state(net, Dict{String,Tuple{Float64,Float64}}(u => (parse(Float64, get(ENV, "PIN", "0")), 1.0) for u in pins), p)
f(res, u) = round(100 * res.node_activities[u], sigdigits = 3)
println("pins: ", [(own[u][1], own[u][2]) for u in pins], "  converged ", r.converged, " (no-pin ", r0.converged, ")")
println("readouts: ", [f(r, u) for u in outs], "  no-pin ", [f(r0, u) for u in outs])
for e in net.edges
    e.edge_type == "context_bridge" || continue
    v, v0 = f(r, e.parent_uuid), f(r0, e.parent_uuid)
    abs(log((v + 1e-6) / (v0 + 1e-6))) > 0.14 || continue
    println("bridge ", own[e.parent_uuid][1], " ", own[e.parent_uuid][2], " -> ", own[e.child_uuid][2], ": ", v, " (no-pin ", v0, ")")
end
