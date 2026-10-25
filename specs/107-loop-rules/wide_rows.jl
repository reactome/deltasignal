using DeltaSignal; const DS = DeltaSignal
for pw in ARGS
    d = "/build/$pw"
    net = redirect_stdout(devnull) do; DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv")); end
    nm(u) = string(net.nodes[u].reactome_id, " ", first(net.nodes[u].display_name, 60))
    m1, _ = DS.apply_material_loops(net; mpools = "on", multistep = "off")
    ws, _ = DS.multistep_carriers(m1, "wide"); ss, _ = DS.multistep_carriers(m1, "strict")
    sset = Set(ss)
    forms = Dict{String, Vector{Tuple{String, Bool}}}()
    for (f, r, c) in (get(ENV,"MODE","wide")=="strict" ? ss : ws); push!(get!(forms, f, Tuple{String,Bool}[]), (r, c)); end
    println("== $pw: wide rows $(length(ws)), strict $(length(ss)), wide-only forms $(length(forms))")
    for f in sort!(collect(keys(forms)), by = nm)
        rs = forms[f]
        println(rpad(nm(f), 72), " kind ", net.node_leaves[f][1], " closed ", rs[1][2], " releases ", length(rs), ": ", join(unique(nm(r)[1:min(end, 45)] for (r, _) in rs[1:min(end, 2)]), " | "))
    end
end
