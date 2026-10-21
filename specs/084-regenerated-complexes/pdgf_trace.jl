using DeltaSignal, CSV, DataFrames, Statistics
const DS = DeltaSignal
d = "/app/pathway_catalog/R-HSA-186797"
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
pins, outs = split.(split(strip(read("/scratch/pdgf.tsv", String)), '\t'), '|')
sid = Dict(r.uuid => r.stable_id for r in eachrow(CSV.read(joinpath(d, "stid_to_uuid_mapping.csv"), DataFrame)))
println("pins: ", [first(split(sid[p], "::")) for p in pins])
for flag in ("off", "on")
    r = withenv(() -> DS.solve_steady_state(net, Dict(String(p) => (0.0, 1.0) for p in pins)), "DS_CARRIER_COMPONENT_TWIN" => flag)
    println("== ", flag, " readout ", round(mean(r.node_activities[u] / 0.01 for u in outs), digits = 3), " component edges ", r.diagnostics["carrier_component_edges"])
    for u in outs
        println("   ", u[1:8], " ", sid[u][1:min(end, 90)], "  ", round(r.node_activities[u] / 0.01, digits = 3), String(u) in pins ? " [pinned]" : "")
    end
end
