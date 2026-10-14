using DeltaSignal; const DS = DeltaSignal
dir = "/catalog/R-HSA-5673001"
net = DS.parse_complete_network(joinpath(dir, "logic_network.csv"), joinpath(dir, "stid_to_uuid_mapping.csv"), nothing)
drugs = DS.drug_uuids(net)
neg = [e for e in net.edges if !e.is_positive && e.parent_uuid in drugs]
println("drug nodes ", length(drugs), "; negative edges from drug nodes ", length(neg),
        " (types ", unique(e.edge_type for e in neg), ")")
net2 = DS.ReactionNetwork(net.nodes, [e for e in net.edges if !(!e.is_positive && e.parent_uuid in drugs)],
                          net.set_mappings, net.cofactor_stids, net.containment, net.drug_stids, net.pools,
                          net.pathogen_stids, net.carrier_loops, net.node_leaves)
cls(v) = v < 0.85 ? 0 : v >= 1.15 ? 2 : 1
p = DS.steady_state_params()
cache = Dict{Tuple{String,String}, Any}()
tally = Dict{String, Int}()
for l in readlines("/s/raf_cases.tsv")
    ax, gene, d, ro, expd, pred0, ui0, pins, outs = split(l, "\t")
    key = (gene, d)
    if !haskey(cache, key)
        v = d == "0" ? 0.0 : 80.0
        obs = Dict{String,Tuple{Float64,Float64}}(String(u) => (v, 1.0) for u in split(pins, "|"))
        r = DS.solve_steady_state(net2, obs, p)
        cache[key] = (r.node_activities, r.converged)
    end
    a, conv = cache[key]
    ui = 100 * sum(a[String(o)] for o in split(outs, "|")) / length(split(outs, "|"))
    pr = cls(ui)
    ok0 = pred0 == expd; ok1 = string(pr) == expd
    k = "$ax $(ok0 ? "R" : "W")->$(ok1 ? "R" : "W")"
    tally[k] = get(tally, k, 0) + 1
    ok0 != ok1 && println(rpad("$ax $gene dir$d $ro", 40), " exp $expd  canon $pred0 ($(round(parse(Float64, ui0), sigdigits=3)))  noDrugReg $pr ($(round(ui, sigdigits=3)))  conv $conv")
end
for k in sort(collect(keys(tally))); println(k, " ", tally[k]); end
