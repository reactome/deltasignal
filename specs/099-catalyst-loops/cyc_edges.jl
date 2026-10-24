using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal
pw, ents = ARGS[1], split(ARGS[2], ',')
d = "/app/pathway_catalog/" * pw
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
sn = DS.solved_carrier_network(net)
nd = CSV.read(joinpath(d, "nodes.csv"), DataFrame); ent = Dict(r.uuid => string(r.diagram_entity_id) for r in eachrow(nd))
sid(u) = get(ent, u, "?")
added = Set((e.parent_uuid, e.child_uuid, e.edge_type) for e in sn.edges) ; orig = Set((e.parent_uuid, e.child_uuid, e.edge_type) for e in net.edges)
for E in ents
    us = [u for u in keys(sn.nodes) if sid(u) == E]
    println("== ", E, " copies ", length(us))
    for u in us
        ins = [(sid(e.parent_uuid), e.edge_type, (e.parent_uuid, e.child_uuid, e.edge_type) in orig ? "" : "SOLVER-ADDED") for e in sn.edges if e.child_uuid == u]
        outs = [(sid(e.child_uuid), e.edge_type) for e in sn.edges if e.parent_uuid == u]
        println("  ", u[1:8], " in: ", ins[1:min(end,6)], "  out: ", length(outs), " ", outs[1:min(end,4)])
    end
end
