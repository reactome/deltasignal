using DeltaSignal; const DS = DeltaSignal
cases = [split(l, "\t") for l in readlines("/s/met_cases.tsv")]
outs = Dict(c[4] => String.(split(c[9], "|")) for c in cases if c[2]=="MET" && c[3]=="0")
met = "93ce1da2-3659-4f60-8060-3ba5ac55c75d"; cbl = "763223ac-59c9-44f7-bc95-17a9b4d1a12d"
hgf = String.(split(first(c for c in cases if c[2]=="HGF")[8], "|"))
for dir in ("/app/pathway_catalog/R-HSA-6806834", "/s/metx")
 net = DS.parse_complete_network(joinpath(dir, "logic_network.csv"), joinpath(dir, "stid_to_uuid_mapping.csv"), nothing)
 for (lab,pins) in (("cblMET",[cbl]), ("bothMET",[cbl, met]), ("HGF",hgf))
  for v in (0.0, 80.0)
   r = DS.solve_steady_state(net, Dict{String,Tuple{Float64,Float64}}(String(u) => (v, 1.0) for u in pins), DS.steady_state_params())
   a = r.node_activities
   m(o) = round(100 * sum(a[u] for u in o) / length(o), sigdigits = 3)
   println(basename(dir), " ", lab, " v=$v ", join(["$(m(o))" for (k,o) in sort(collect(outs))], " "), " conv=", r.converged)
  end
 end
end
