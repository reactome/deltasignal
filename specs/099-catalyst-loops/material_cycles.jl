# specs/099: the shortest material cycle of each unhandled material-closed positive
# loop (specs/097's definition), with the leaves of every entity on it, for
# classify_cycles.py. Usage (dev container): julia material_cycles.jl BUILD OUT.tsv
using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal
include(joinpath(@__DIR__, "..", "097-loop-audit", "sccs.jl"))
build, out_path = ARGS[1], ARGS[2]
const UB = Set(["R-HSA-68524", "R-HSA-113595", "R-HSA-9660007"])
rows = DataFrame(pathway = String[], comp_size = Int[], cycle_len = Int[], entities = String[], leaves = String[])
function shortest_cycle(nodes, adj)
    best = String[]
    for s in nodes[1:min(end, 40)]                     # a few starts suffice for the shortest
        prev = Dict{String, String}(); q = [s]; seen = Set([s]); found = ""
        while !isempty(q) && isempty(found)
            v = popfirst!(q)
            for w in get(adj, v, String[])
                w == s && (found = v; break)
                w in seen || (push!(seen, w); prev[w] = v; push!(q, w))
            end
        end
        isempty(found) && continue
        cyc = [found]; while cyc[end] != s; push!(cyc, prev[cyc[end]]); end
        (isempty(best) || length(cyc) < length(best)) && (best = reverse(cyc))
    end
    best
end
for pw in sort(filter(p -> startswith(p, "R-HSA-"), readdir(build)))
    d = joinpath(build, pw); isfile(joinpath(d, "logic_network.csv")) || continue
    net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
    sn = DS.solved_carrier_network(net)
    nd = CSV.read(joinpath(d, "nodes.csv"), DataFrame; types = String)
    kind = Dict(r.uuid => r.node_kind for r in eachrow(nd))
    leaves = Dict(r.uuid => coalesce(r.member_leaves, "") for r in eachrow(nd))
    sid(u) = (n = get(sn.nodes, u, nothing); n === nothing ? "" : something(n.reactome_id, ""))
    pool_of = Dict{String, Int}()
    if net.pools !== nothing
        for (k, p) in enumerate(net.pools)
            for u in vcat(p.forms, p.intermediates); pool_of[u] = k; end
            for q in p.paths, (a, b, rxs) in q.steps, rx in rxs; pool_of[rx] = k; end
        end
    end
    self_inh = first(DS.self_inhibitor_setup(sn)); self_act = first(DS.self_activator_setup(sn))
    adj = Dict{String, Vector{String}}(); mat = Dict{String, Vector{String}}(); neg = Set{Tuple{String, String}}()
    for e in sn.edges
        a, b = e.parent_uuid, e.child_uuid
        (haskey(pool_of, a) && get(pool_of, b, 0) == pool_of[a]) && continue
        ((!e.is_positive && haskey(self_inh, (a, b))) || (e.is_positive && haskey(self_act, (a, b)))) && continue
        push!(get!(adj, a, String[]), b)
        e.is_positive || push!(neg, (a, b))
        e.edge_type in ("catalyst", "regulator") || push!(get!(mat, a, String[]), b)
    end
    for c in sccs(sort(collect(keys(sn.nodes))), adj)
        cs = Set(c)
        any(((a, b),) -> (a, b) in neg, [(a, b) for a in c for b in get(adj, a, String[]) if b in cs]) && continue
        madj = Dict(a => [b for b in get(mat, a, String[]) if b in cs] for a in c)
        mcs = sccs(sort(c), madj)
        isempty(mcs) && continue                                  # catalyst-closed: not here
        any(u -> sid(u) in UB, reduce(vcat, mcs)) && continue     # ubiquitin-closed: not here
        for mc in mcs
            mset = Set(mc)
            sub = Dict(a => [b for b in madj[a] if b in mset] for a in mc)
            cyc = shortest_cycle(sort(mc), sub)
            ents = [u for u in cyc if get(kind, u, "") != "reaction"]
            push!(rows, (pw, length(c), length(cyc), join([sid(u) for u in ents], "|"),
                         join([get(leaves, u, "") for u in ents], ";")))
        end
    end
end
CSV.write(out_path, rows; delim = '\t')
println("material cycles: ", nrow(rows), " in ", length(unique(rows.pathway)), " pathways")
