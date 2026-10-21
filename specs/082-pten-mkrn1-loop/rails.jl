using DeltaSignal, CSV, DataFrames, Statistics
const DS = DeltaSignal
const PW = Dict("Signaling_by_Activin" => "R-HSA-1502540", "Signaling_by_EGFR" => "R-HSA-177929",
                "RAF_MAP_kinase_cascade" => "R-HSA-5673001", "Cell_Cycle_Checkpoints" => "R-HSA-69620")
for line in eachline("/scratch/rails.tsv")
    pw, gene, ko, _, pins, outs = split(line, '\t')
    d = "/app/pathway_catalog/" * PW[pw]
    net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
    P = split(pins, '|'); O = split(outs, '|')
    nd = CSV.read(joinpath(d, "nodes.csv"), DataFrame); ent = Dict(r.uuid => string(r.diagram_entity_id) for r in eachrow(nd))
    kind = Dict(r.uuid => r.node_kind for r in eachrow(nd))
    vals = String[]
    res = Dict{Float64, Any}()
    for v in (0.0, 0.5, 0.9, 0.99, 1.0, 1.01, 1.1, 2.0, 80.0)
        r = DS.solve_steady_state(net, Dict(String(p) => (v, 1.0) for p in P)); res[v] = r
        push!(vals, string(v, "→", round(mean(r.node_activities[u] / 0.01 for u in O), digits = 3)))
    end
    println("\n=== ", pw, " ", gene, " readout ", ko, ": ", join(vals, "  "))
    a, b = res[1.0], res[1.01]
    rows = [(u, b.node_activities[u] / 0.01) for u in keys(net.nodes) if get(kind, u, "") != "reaction"]
    sort!(rows, by = x -> -abs(log(max(x[2], 1e-6))))
    seen = Set{String}()
    for (u, y) in rows
        ent[u] in seen && continue; push!(seen, ent[u]); length(seen) > 10 && break
        println("   ", rpad(ent[u], 15), " at 1.01×: ", round(y, digits = 4))
    end
end
