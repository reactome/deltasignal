# specs/061 R1 trace: TP53 pathway, MDM2 knockout (pins as the rtwin58g case
# file lists them), under a given DS_DEPLETION_OWN_PRODUCT; prints the source
# copies of each entity-own pair and the depleters.
using DeltaSignal; const DS = DeltaSignal
dir = ENV["PW"]
net = DS.parse_complete_network(joinpath(dir, "logic_network.csv"), joinpath(dir, "stid_to_uuid_mapping.csv"), nothing)
pins = split(ENV["PINS"], "|")
r = DS.solve_steady_state(net, Dict{String,Tuple{Float64,Float64}}(String(u) => (0.0, 1.0) for u in pins), DS.steady_state_params())
println("mode ", get(ENV, "DS_DEPLETION_OWN_PRODUCT", "full"), " conv ", r.converged, " iters ", r.iterations)
src = DS.own_source_depletion_sources(net)
for ((d, t), ss) in sort(collect(src))
    println("pair D=", d[1:8], " ", net.nodes[d].reactome_id, " fold ", round(100r.node_activities[d], sigdigits=4),
            "  T=", t[1:8], " ", net.nodes[t].reactome_id, " fold ", round(100r.node_activities[t], sigdigits=4),
            "  sources ", join([s[1:8] * "=" * string(round(100r.node_activities[s], sigdigits=4)) for s in ss], ","))
end
for u in split(get(ENV, "SHOW", ""), "|")
    isempty(u) || println("show ", u, " ", net.nodes[String(u)].reactome_id, " ", 100r.node_activities[String(u)])
end
