# For each discordant case: solve under the defaults at several iteration caps;
# report the predicted class per cap and convergence. Reads /probe/disc.tsv.
using DeltaSignal
cls(v) = v < 0.85 ? "0" : (v > 1.15 ? "2" : "1")
rows = [split(l, '\t') for l in eachline("/probe/disc.tsv")]
nets = Dict{String,Any}()
caps = [500, 1000, 2000, 4000]
open("/probe/capcheck.tsv", "w") do io
  for r in rows
    ax, split_, pw, pid, gene, dir, ko, expd, cpred, apred, apui, pins, outs = r
    net = get!(nets, pid) do
      d = "/app/pathway_catalog/" * pid
      DeltaSignal.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"), nothing; boundary_path=:sibling)
    end
    val = dir == "0" ? 0.0 : 80.0
    obs = Dict{String,Tuple{Float64,Float64}}(String(p) => (val, 1.0) for p in split(pins, '|') if !isempty(p))
    ou = [String(u) for u in split(outs, '|') if !isempty(u)]
    res = String[]
    for mi in caps
      ENV["DS_MAX_ITERS"] = string(mi)
      s = DeltaSignal.solve_steady_state(net, obs)
      v = sum(get(s.node_activities, u, 0.01) for u in ou) / length(ou) * 100
      push!(res, string(cls(v), ":", round(v, digits = 3), ":", s.converged ? "c" : "n"))
    end
    println(io, join([ax, split_, pw, gene, dir, ko, expd, cpred, apred, apui, res...], '\t'))
    flush(io)
  end
end
