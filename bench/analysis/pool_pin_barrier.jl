#!/usr/bin/env julia
# specs/052 R2 (PTK6): a pin is a reachability barrier in the cycle-balance
# pool's supply rule (specs/039), so pinning a node on a pool's own recycling
# path turns that recycling step into "external supply".
#
# Case: PTK6 signalling (R-HSA-8848021), DOK1 knockout. Every pin below is set
# at the value the unpinned solve already gives that node, so the pin changes
# nothing except that the node is pinned.
#
# HOW TO RUN (dev container, canonical build mounted at /catalog):
#   docker compose -f docker-compose.dev.yml run --rm --no-deps \
#     -v ~/deltasignal-catalogs/builds/20261012-1401_01206ca:/catalog:ro \
#     julia-api julia --project=/app /app/bench/analysis/pool_pin_barrier.jl
#
# Expected (solver at the commit that added this file, code defaults):
#   no extra pin                         PTK6 state 1.0       nodes moved 0
#   pin DOK1:PTK6 at its own value       PTK6 state 0.888889  nodes moved 269
#   pin R-HSA-8848776 at its own value   PTK6 state 0.888892  nodes moved 268
# and, sweeping the DOK1:PTK6 pin from 1e-3 to exactly 0, 0.888889 every time.
using DeltaSignal
const DS = DeltaSignal

const DIR = get(ENV, "PTK6_DIR", "/catalog/R-HSA-8848021")
net = DS.parse_complete_network(joinpath(DIR, "logic_network.csv"),
                                joinpath(DIR, "stid_to_uuid_mapping.csv"), nothing)
node(prefix) = only(u for u in keys(net.nodes) if startswith(u, prefix))
const DOK1 = node("c53d802c")        # free DOK1, the knockout pin
const CX = node("dacc03ed")          # DOK1:PTK6, R-HSA-8848773
const RX = node("b69fd0b2")          # R-HSA-8848776: consumes DOK1:PTK6, releases PTK6
const STATE = node("c38718c8")       # PTK6 form R-HSA-8848007, a state of pool1
const POOL = "pool1"

p = DS.steady_state_params()
ko = Dict{String, Tuple{Float64, Float64}}(DOK1 => (0.0, 1.0))
base = DS.solve_steady_state(net, ko, p).node_activities

function report(label, extra::Dict{String, Float64})
    obs = copy(ko)
    for (u, v) in extra
        obs[u] = (v, 1.0)
    end
    a = DS.solve_steady_state(net, obs, p).node_activities
    moved = count(u -> abs(a[u] - base[u]) > 1e-9, keys(base))
    println(rpad(label, 44), " PTK6 state ", round(a[STATE] * 100, sigdigits = 6), "  nodes moved ", moved)
end

report("no extra pin", Dict{String, Float64}())
report("pin DOK1:PTK6 at its own value", Dict(CX => base[CX] * 100))
report("pin R-HSA-8848776 at its own value", Dict(RX => base[RX] * 100))
for v in (1e-3, 1e-6, 1e-10, 1e-14, 0.0)
    report("pin DOK1:PTK6 at $v", Dict(CX => v))
end

# Which non-transition producers of the pool's states the pool reaches, with
# and without DOK1:PTK6 as a barrier (positive edges inside the component;
# an approximation of the solver's pool_out, enough to name the producer).
ids = collect(keys(net.nodes)); ix = Dict(u => i for (i, u) in enumerate(ids))
adj = [Set{Int}() for _ in ids]
for e in net.edges
    haskey(ix, e.parent_uuid) && haskey(ix, e.child_uuid) && push!(adj[ix[e.parent_uuid]], ix[e.child_uuid])
end
comp, _ = DS.tarjan_scc_components(adj)
rows(f) = [split(l, ",") for l in eachline(joinpath(DIR, f)) if startswith(l, POOL * ",")]
forms = [r[2] for r in rows("pools.csv") if r[4] == "state"]
steps = Set(r[6] for r in rows("pool_transitions.csv"))
prods = sort(unique(e.parent_uuid for e in net.edges
                    if e.child_uuid in forms && e.is_positive && !(e.parent_uuid in steps)))
function reach(barrier)
    s = Set(forms); q = copy(forms)
    while !isempty(q)
        v = pop!(q)
        for e in net.edges
            (e.parent_uuid == v && e.is_positive) || continue
            w = e.child_uuid
            (w in s || w in barrier || comp[ix[w]] != comp[ix[STATE]]) && continue
            push!(s, w); push!(q, w)
        end
    end
    s
end
r0, r1 = reach(Set([DOK1])), reach(Set([DOK1, CX]))
println("\nsupply producers of $POOL (not reached from its own states):")
println("  unpinned:        ", count(u -> !(u in r0), prods), " of ", length(prods))
println("  DOK1:PTK6 pinned: ", count(u -> !(u in r1), prods), " of ", length(prods))
for u in prods
    (u in r0) == (u in r1) && continue
    println("  becomes supply when pinned: ", net.nodes[u].reactome_id, " ", u[1:8],
            ", reads ", round(base[u] * 100, sigdigits = 4))
end
