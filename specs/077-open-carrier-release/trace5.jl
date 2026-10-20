using DeltaSignal, CSV, DataFrames, Statistics
const DS = DeltaSignal
const PW = Dict("S_Phase" => "R-HSA-69242", "Signaling_by_ERBB2" => "R-HSA-1227986")
for line in eachline("ARGS[1]")
    f = split(line, '\t'); pw, gene, dir, ro = f[1:4]
    d = "/app/pathway_catalog/" * PW[pw]
    net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
    pins = split(f[8], '|'); outs = split(f[9], '|')
    obs = Dict(String(p) => (dir == "0" ? 0.0 : 80.0, 1.0) for p in pins)
    r = Dict(o => withenv(() -> DS.solve_steady_state(net, obs), "DS_CARRIER_LOOP_OPEN" => o) for o in ("off", "fed_single"))
    v(o, u) = r[o].node_activities[u] / 0.01
    nd = CSV.read(joinpath(d, "nodes.csv"), DataFrame); ent = Dict(x.uuid => string(x.diagram_entity_id) for x in eachrow(nd))
    println("\n=== $pw $gene $dir -> $ro: readout off ", round(mean(v.("off", outs)), digits=3), " fed_single ", round(mean(v.("fed_single", outs)), digits=3))
    sn = withenv(() -> DS.solved_carrier_network(net), "DS_CARRIER_LOOP_OPEN" => "fed_single")
    kept = Set((e.parent_uuid, e.child_uuid) for e in sn.edges)
    seen = Set{String}()
    for e in net.edges
        (e.parent_uuid, e.child_uuid) in kept && continue
        k = ent[e.parent_uuid] * "->" * ent[e.child_uuid]
        k in seen && continue; push!(seen, k)
        c = e.child_uuid
        abs(log(max(v("off", c),1e-6)) - log(max(v("fed_single", c),1e-6))) > 0.05 &&
            println("  dropped ", k, "  child off ", round(v("off", c), digits=3), " fed_single ", round(v("fed_single", c), digits=3))
    end
    rows = [(u, v("off", u), v("fed_single", u)) for u in keys(net.nodes) if haskey(r["off"].node_activities, u)]
    filter!(x -> abs(log(max(x[2],1e-6)) - log(max(x[3],1e-6))) > 0.1, rows)
    sort!(rows, by = x -> -abs(log(max(x[2],1e-6)) - log(max(x[3],1e-6))))
    es = Set{String}()
    for (u, a, b) in rows
        ent[u] in es && continue; push!(es, ent[u]); length(es) > 12 && break
        println("    ", ent[u], " ", u[1:8], " off ", round(a, digits=3), " fed_single ", round(b, digits=3), u in outs ? "  [readout]" : "")
    end
end
