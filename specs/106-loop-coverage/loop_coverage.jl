# specs/106: loop coverage. Every loop in the network AS GENERATED, the solver
# rules that handle edges inside it, and what is left of it once they have.
#
# Mechanisms (each as the default solver applies it):
#   recycle - recycle edges (specs/081), solved by the pool balance;
#   carrier - carrier-loop edges the carrier rules drop (specs/049, 077, 088);
#   pool    - an edge between two nodes of one pool (specs/039, 081, 086);
#   self    - an inhibitor or activator slot built from its own reaction's
#             input (specs/022, 040, 083, 095).
# A loop's remnants are the cyclic components left inside it on the solved
# network (carrier edges dropped, twins and supply edges added, recycle edges
# stripped) once pool and self edges are removed, labelled as in specs/097.
#
# Usage (dev container): julia loop_coverage.jl BUILD_DIR OUT_PREFIX
using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal

build, prefix = ARGS[1], ARGS[2]
const UB = Set(["R-HSA-68524", "R-HSA-113595", "R-HSA-9660007"])   # ubiquitin, as DS.UBIQUITIN_STIDS

function sccs(nodes::Vector{String}, adj::Dict{String, Vector{String}})
    # iterative Tarjan
    index = Dict{String, Int}(); low = Dict{String, Int}(); onstack = Set{String}()
    stack = String[]; out = Vector{Vector{String}}(); i = 0
    for root in nodes
        haskey(index, root) && continue
        work = [(root, 1)]
        while !isempty(work)
            v, k = pop!(work)
            if k == 1
                i += 1; index[v] = i; low[v] = i; push!(stack, v); push!(onstack, v)
            end
            nbrs = get(adj, v, String[])
            recurse = false
            for j in k:length(nbrs)
                w = nbrs[j]
                if !haskey(index, w)
                    push!(work, (v, j + 1)); push!(work, (w, 1)); recurse = true; break
                elseif w in onstack
                    low[v] = min(low[v], index[w])
                end
            end
            recurse && continue
            if low[v] == index[v]
                comp = String[]
                while true
                    w = pop!(stack); delete!(onstack, w); push!(comp, w)
                    w == v && break
                end
                length(comp) > 1 && push!(out, comp)
            end
            if !isempty(work)
                u = work[end][1]
                low[u] = min(low[u], low[v])
            end
        end
    end
    out
end

function reach(adj, from)
    seen = Set(from); st = copy(from)
    while !isempty(st)
        v = pop!(st)
        for w in get(adj, v, String[]); w in seen || (push!(seen, w); push!(st, w)); end
    end
    seen
end


loop_rows = DataFrame(pathway = String[], size = Int[], touched = String[], status = String[],
                      remnants = String[], remnant_nodes = Int[], example = String[])
pw_rows = DataFrame(pathway = String[], nodes = Int[], loops = Int[], loop_nodes = Int[],
                    recycle_edges = Int[], carrier_edges = Int[], pool_edges = Int[], self_edges = Int[],
                    created_remnants = Int[])

for pw in sort(filter(p -> startswith(p, "R-HSA-"), readdir(build)))
    d = joinpath(build, pw)
    isfile(joinpath(d, "logic_network.csv")) || continue
    net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
    sn = DS.solved_carrier_network(net)
    sid(u) = (n = get(net.nodes, u, nothing); n === nothing ? "" : something(n.reactome_id, ""))
    # handled edges by mechanism
    recycle = Set(last(DS.split_recycle_edges(net)))
    key(e) = (e.parent_uuid, e.child_uuid, e.edge_type, e.is_positive)
    kept = Set(key(e) for e in sn.edges)
    carrier = Set((e.parent_uuid, e.child_uuid) for e in net.edges
                  if e.edge_type != DS.RECYCLE_EDGE && !(key(e) in kept))
    pool_of = Dict{String, Int}()
    if net.pools !== nothing
        for (k, p) in enumerate(net.pools)
            for u in vcat(p.forms, p.intermediates); pool_of[u] = k; end
            for q in p.paths, (a, b, rxs) in q.steps, rx in rxs; pool_of[rx] = k; end
        end
    end
    self_inh = first(DS.self_inhibitor_setup(sn))
    self_act = first(DS.self_activator_setup(sn))
    pool = Set{Tuple{String, String}}(); slf = Set{Tuple{String, String}}()
    adj = Dict{String, Vector{String}}(); neg = Set{Tuple{String, String}}()
    etype = Dict{Tuple{String, String}, String}()
    for e in sn.edges
        a, b = e.parent_uuid, e.child_uuid
        etype[(a, b)] = e.edge_type
        e.is_positive || push!(neg, (a, b))
        if haskey(pool_of, a) && get(pool_of, b, 0) == pool_of[a]
            push!(pool, (a, b)); continue
        end
        if (!e.is_positive && haskey(self_inh, (a, b))) || (e.is_positive && haskey(self_act, (a, b)))
            push!(slf, (a, b)); continue
        end
        push!(get!(adj, a, String[]), b)
    end
    raw_adj = Dict{String, Vector{String}}()
    for e in net.edges; push!(get!(raw_adj, e.parent_uuid, String[]), e.child_uuid); end
    nodes = sort(unique(vcat(collect(keys(net.nodes)), collect(keys(sn.nodes)))))
    loops = sccs(nodes, raw_adj)
    loop_of = Dict(u => k for (k, c) in enumerate(loops) for u in c)
    # remnants, labelled as in specs/097
    rem = Dict{Int, Vector{Tuple{String, Int}}}(); created = 0
    for c in sccs(nodes, adj)
        cs = Set(c)
        inner = [(a, b) for a in c for b in get(adj, a, String[]) if b in cs]
        nneg = count(in(neg), inner)
        selfdep = any(((a, b),) -> etype[(a, b)] == "depletion" &&
                      sid(b) in get(sn.containment, sid(a), Set{String}()), inner)
        lab = nneg == 0 ? "positive" : selfdep ? "self_depleting" : "negative"
        if lab == "positive"
            madj = Dict{String, Vector{String}}()
            for (a, b) in inner
                etype[(a, b)] in ("catalyst", "regulator") && continue
                push!(get!(madj, a, String[]), b)
            end
            mc = sccs(sort(c), madj)
            lab *= ":" * (isempty(mc) ? "catalyst" : any(u -> sid(u) in UB, reduce(vcat, mc)) ? "ubiquitin" : "material")
        end
        owners = unique([get(loop_of, u, 0) for u in c])
        if length(owners) == 1 && owners[1] > 0
            push!(get!(rem, owners[1], Tuple{String, Int}[]), (lab, length(c)))
        else
            created += 1   # closed by an edge the solver adds, or spanning two generated loops
        end
    end
    inside(k, s) = count(((a, b),) -> get(loop_of, a, 0) == k && get(loop_of, b, 0) == k, s)
    for (k, c) in enumerate(loops)
        t = String[]
        inside(k, recycle) > 0 && push!(t, "recycle")
        inside(k, carrier) > 0 && push!(t, "carrier")
        inside(k, pool) > 0 && push!(t, "pool")
        inside(k, slf) > 0 && push!(t, "self")
        r = get(rem, k, Tuple{String, Int}[])
        status = isempty(r) ? (isempty(t) ? "acyclic_after_solver" : "dissolved") :
                 isempty(t) ? "untouched" : "partly_handled"
        push!(loop_rows, (pw, length(c), join(t, "|"), status, join(sort(unique(first.(r))), "+"),
                          sum(last, r; init = 0), sid(first(sort(c)))))
    end
    push!(pw_rows, (pw, length(nodes), length(loops), sum(length, loops; init = 0),
                    length(recycle), length(carrier), length(pool), length(slf), created))
end
CSV.write(prefix * "_pathways.tsv", pw_rows; delim = '\t')
CSV.write(prefix * "_loops.tsv", loop_rows; delim = '\t')
println("pathways ", nrow(pw_rows), "; with loops ", count(>(0), pw_rows.loops), "; loops ", nrow(loop_rows),
        "; remnants closed by solver-added edges ", sum(pw_rows.created_remnants))
for s in ("dissolved", "partly_handled", "untouched", "acyclic_after_solver")
    x = loop_rows[loop_rows.status .== s, :]
    println("  $s: ", nrow(x), " loops, ", sum(x.size; init = 0), " nodes, ", length(unique(x.pathway)), " pathways")
end
for m in ("recycle", "carrier", "pool", "self")
    x = loop_rows[occursin.(m, loop_rows.touched), :]
    println("  touched by $m: ", nrow(x), " loops, ", length(unique(x.pathway)), " pathways")
end
