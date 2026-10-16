using DeltaSignal
d = "/app/pathway_catalog/R-HSA-3700989"
net = DeltaSignal.parse_complete_network(joinpath(d,"logic_network.csv"), joinpath(d,"stid_to_uuid_mapping.csv"), nothing; boundary_path=:sibling)
for pinset in ("adopt",)
  l = readlines("/probe/tp53_$(pinset).txt")
  pins = [String(p) for p in split(l[1], '|') if !isempty(p)]; outs = [String(u) for u in split(l[2], '|') if !isempty(u)]
  obs = Dict{String,Tuple{Float64,Float64}}(p => (0.0, 1.0) for p in pins)
  for (own, twin) in (("suppress_only","off"), ("suppress_only_entity","off"))
    ENV["DS_DEPLETION_OWN_PRODUCT"] = own; ENV["DS_RECYCLED_TWIN"] = twin
    vals = String[]
    for mi in (500, 2000)
      ENV["DS_MAX_ITERS"] = string(mi)
      r = DeltaSignal.solve_steady_state(net, obs)
      v = sum(get(r.node_activities, u, 0.01) for u in outs)/length(outs)*100
      push!(vals, string(round(v, digits=3), r.converged ? "c" : "n", "(", round(r.final_residual, sigdigits=2), ")"))
    end
    println(pinset, " pins=", length(pins), " own=", own, " twin=", twin, " ", join(vals, "  "))
  end
end
