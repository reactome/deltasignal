using DeltaSignal, CSV, DataFrames, Statistics
const DS = DeltaSignal
d = "/app/pathway_catalog/R-HSA-453279"
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
nd = CSV.read(joinpath(d, "nodes.csv"), DataFrame)
kind = Dict(r.uuid => r.node_kind for r in eachrow(nd)); ent = Dict(r.uuid => string(r.diagram_entity_id) for r in eachrow(nd))
produced = Set(e.child_uuid for e in net.edges if e.is_positive && e.edge_type == "output")
leaf = Set(e.parent_uuid for e in net.edges if e.edge_type == "assembly")
hasout = Set(e.parent_uuid for e in net.edges)
byent = Dict{String,Vector{String}}()
for (u, _) in net.nodes; push!(get!(byent, get(ent, u, ""), String[]), u); end
add = DS.LogicNetworkEdge[]
for (en, us) in byent
    isempty(en) && continue
    ps = [u for u in us if u in produced && get(kind, u, "") != "dissociation_sink" &&
          (get(ENV, "DEADEND", "0") == "0" || !(u in hasout))]
    bs = [u for u in us if u in leaf]
    for p in ps, b in bs
        p == b && continue
        push!(add, DS.LogicNetworkEdge(p, b, false, true, 1.0, "produced_supply"))
    end
end
println("added edges: ", length(add), "  ", sort(unique(ent[e.parent_uuid] for e in add)))
net2 = DS.ReactionNetwork(net.nodes, vcat(net.edges, add), net.set_mappings, net.cofactor_stids, net.containment,
                          net.drug_stids, net.pools, net.pathogen_stids, net.carrier_loops, net.node_leaves)
base2 = DS.solve_steady_state(net2, Dict{String,Tuple{Float64,Float64}}())
off = count(u -> abs(base2.node_activities[u] / 0.01 - 1) > 0.01, collect(keys(base2.node_activities)))
println("unperturbed probe network: nodes off baseline by >1%: ", off, " converged ", base2.converged)
cls(x) = x < 0.85 ? 0 : (x >= 1.15 ? 2 : 1)
rows = [split(l, '\t') for l in eachline(ARGS[1])]
memo = Dict()
tally = Dict{String,Int}()
for f in rows
    gene, dir, ko, expd, pred0, pins, outs = f
    key = (gene, dir)
    if !haskey(memo, key)
        obs = Dict(String(p) => (dir == "0" ? 0.0 : 80.0, 1.0) for p in split(pins, '|'))
        memo[key] = (DS.solve_steady_state(net, obs), DS.solve_steady_state(net2, obs))
    end
    a, b = memo[key]
    o = split(outs, '|')
    va = mean(a.node_activities[u] / 0.01 for u in o); vb = mean(b.node_activities[u] / 0.01 for u in o)
    ca, cb = cls(va), cls(vb)
    string(ca) != pred0 && (tally["baseline mismatch"] = get(tally, "baseline mismatch", 0) + 1)
    k = string(ca == parse(Int, expd) ? "right" : "wrong", " -> ", cb == parse(Int, expd) ? "right" : "wrong")
    tally[k] = get(tally, k, 0) + 1
    if ca != cb && length(ARGS) > 1
        println(gene, " ", dir, " ", ko, " exp ", expd, "  base ", round(va, digits=3), "  probe ", round(vb, digits=3), " conv ", b.converged)
    end
end
println(tally)
