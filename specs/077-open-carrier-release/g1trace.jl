using DeltaSignal, CSV, DataFrames, Statistics
const DS = DeltaSignal
d = "/app/pathway_catalog/R-HSA-453279"
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
pins = ["9dee8c84-63ff-4903-a7e9-d1335ccae245", "b77963c4-f2be-4fe1-919c-768a6de5c50f"]
outs = split(read("ARGS[1]", String) |> strip, "|")
obs = Dict(p => (80.0, 1.0) for p in pins)
function go(open)
    withenv("DS_CARRIER_LOOP_OPEN" => open) do
        DS.solve_steady_state(net, obs)
    end
end
a = go("off"); b = go("fed_single")
fa(u) = a.node_activities[u] / 0.01; fb(u) = b.node_activities[u] / 0.01
println("readout mean off ", mean(fa.(outs)), " fed_single ", mean(fb.(outs)), " conv ", a.converged, " ", b.converged)
nd = CSV.read(joinpath(d, "nodes.csv"), DataFrame)
ent = Dict(r.uuid => string(r.diagram_entity_id) for r in eachrow(nd))
rows = [(u, fa(u), fb(u)) for u in keys(a.node_activities) if abs(log(max(fa(u),1e-9)) - log(max(fb(u),1e-9))) > 0.5]
sort!(rows, by = r -> -abs(log(max(r[2],1e-9)) - log(max(r[3],1e-9))))
println("nodes differing (", length(rows), ")")
for (u, x, y) in rows[1:min(end, 60)]
    println(rpad(get(ent, u, "?"), 16), " ", u[1:8], "  off ", round(x, digits=3), "  fed_single ", round(y, digits=3))
end
dr = Set((e.parent_uuid, e.child_uuid) for e in net.edges) 
open_sn = withenv(() -> DS.solved_carrier_network(net), "DS_CARRIER_LOOP_OPEN" => "fed_single")
kept = Set((e.parent_uuid, e.child_uuid) for e in open_sn.edges)
println("dropped edges:")
for (p, c) in setdiff(dr, kept)
    println("  ", get(ent, p, "?"), " ", p[1:8], " -> ", get(ent, c, "?"), " ", c[1:8], "  off ", round(fa(c),digits=3), " fed_single ", round(fb(c),digits=3))
end
