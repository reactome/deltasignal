using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal
d = "/app/pathway_catalog/R-HSA-453279"
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
nd = CSV.read(joinpath(d, "nodes.csv"), DataFrame); ent = Dict(r.uuid => string(r.diagram_entity_id) for r in eachrow(nd))
pin = first(u for u in keys(net.nodes) if startswith(u, "3cd000b2"))
pl = only(p for p in net.pools if p.id == "recycle2")
for m in ("off", "on"), v in (80.0, 0.0)
    r = withenv("DS_POOL_SUPPLY_AND" => m) do
        DS.solve_steady_state(net, Dict(pin => (v, 1.0)))
    end
    f(u) = round(r.node_activities[u] / 0.01, digits = 3)
    println("== $m RBL2=$v converged=$(r.converged) resid=$(r.final_residual)")
    println("  pool forms: ", [(ent[u][1:min(end,40)], f(u)) for u in pl.forms])
    println("  CCNA1 68891: ", sort(unique([f(u) for u in keys(net.nodes) if ent[u] == "R-HSA-68891"])))
    # nodes read from pool forms
    kids = Set(e.child_uuid for e in net.edges if e.parent_uuid in pl.forms && !(e.child_uuid in pl.forms))
    println("  consumers: ", sort([(ent[u][1:min(end,30)], f(u)) for u in kids])[1:min(end,12)])
end
