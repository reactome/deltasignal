# specs/056 check: a pruned context bundle gives the pathway's own nodes the
# same values as the unpruned one. Pins each of up to 12 own roots to 0 and 80.
using DeltaSignal; const DS = DeltaSignal
load(d) = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"), nothing)
p = DS.steady_state_params()
for P in ARGS
    a, b = load("/ctxp/$P"), load("/ctxu/$P")
    own = [l for l in readlines("/ctxp/$P/context_own.csv")[2:end]]
    indeg = Dict{String,Int}()
    for e in a.edges; indeg[e.child_uuid] = get(indeg, e.child_uuid, 0) + 1; end
    roots = sort([u for u in own if get(indeg, u, 0) == 0 && haskey(a.nodes, u)])[1:min(12, end)]
    worst = 0.0; nconv = 0
    for u in roots, v in (0.0, 80.0)
        obs = Dict{String,Tuple{Float64,Float64}}(u => (v, 1.0))
        ra, rb = DS.solve_steady_state(a, obs, p), DS.solve_steady_state(b, obs, p)
        nconv += ra.converged == rb.converged
        worst = max(worst, maximum(abs(ra.node_activities[w] - rb.node_activities[w]) for w in own if haskey(ra.node_activities, w)))
    end
    println(P, ": nodes ", length(a.nodes), " vs ", length(b.nodes), "; ", 2length(roots), " solves; max |diff| on own nodes ", worst, "; converged flag equal ", nconv)
end
