using DeltaSignal, CSV, DataFrames, Statistics
const DS = DeltaSignal
d = "/app/pathway_catalog/R-HSA-1257604"
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
pins, outs = split.(split(strip(read("/scratch/pten.tsv", String)), '\t'), '|')
nd = CSV.read(joinpath(d, "nodes.csv"), DataFrame); ent = Dict(r.uuid => string(r.diagram_entity_id) for r in eachrow(nd))
println("PTEN pins: ", [ent[p] for p in pins])
for v in (0.0, 0.5, 1.0, 1.1, 2.0, 10.0, 80.0)
    r = DS.solve_steady_state(net, Dict(String(p) => (v, 1.0) for p in pins))
    println(rpad(v, 6), " readout ", round(mean(r.node_activities[u] / 0.01 for u in outs), digits = 3), "  converged ", r.converged)
end
r0 = DS.solve_steady_state(net, Dict{String, Tuple{Float64, Float64}}())
println("no pins: readout ", round(mean(r0.node_activities[u] / 0.01 for u in outs), digits = 3))
