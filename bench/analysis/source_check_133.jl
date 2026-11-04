# specs/133 amendment 2: does a source pathway inside a one-hop context bundle solve exactly as it does alone?
# For every bundle pathway and every source pathway Q in it, pins each of up to K of Q's starting copies (sorted by
# uuid; an in-edge-free node with an out-edge) at 5x, one at a time, and compares every node of Q in the bundle with Q
# solved alone from the base build (same uuids). Under DS_CONTEXT_BRIDGE_ORDER=last, as the arms run.
# Usage (in the dev container): julia --project=/app source_check_133.jl BUNDLE_DIR BASE_DIR [K]
using DeltaSignal
const D = DeltaSignal
ENV["DS_CONTEXT_BRIDGE_ORDER"] = "last"
bundle, base = ARGS[1], ARGS[2]
K = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 3
load(d) = D.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
alone = Dict{String, Any}()
worst_all = 0.0
for p in sort(readdir(bundle))
    d = joinpath(bundle, p)
    isfile(joinpath(d, "context_origin.csv")) || continue
    net = load(d)
    origin = net.node_origin
    qs = sort(unique(q for q in values(origin) if q != p))
    isempty(qs) && continue
    for q in qs
        qnet = get!(() -> load(joinpath(base, q)), alone, q)
        hasin = Set(e.child_uuid for e in qnet.edges)
        hasout = Set(e.parent_uuid for e in qnet.edges)
        starts = sort([u for u in keys(qnet.nodes) if !(u in hasin) && u in hasout && get(origin, u, "") == q])
        worst = 0.0
        for u in starts[1:min(K, length(starts))]
            obs = Dict(u => (5.0, 1.0))
            a = D.solve_steady_state(net, obs).node_activities
            b = D.solve_steady_state(qnet, obs).node_activities
            for (v, x) in b
                haskey(a, v) || continue
                worst = max(worst, abs(log((a[v] + 1e-12) / (x + 1e-12))))
            end
        end
        global worst_all = max(worst_all, worst)
        println(p, "\t", q, "\tpins=", min(K, length(starts)), "\tmax_abs_log_ratio=", round(worst, sigdigits = 3))
    end
end
println("worst over all: ", worst_all)
