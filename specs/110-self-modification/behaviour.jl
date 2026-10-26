# specs/110 behaviour checks (B2-B5), specs/107's harness with a new selection.
#
#   julia behaviour.jl BUILD select OUT.tsv
#       per pathway, the nodes of the cyclic components that contain a node which is
#       both an input and a catalyst of a reaction node in the same component, on the
#       network as the solver sees it with every DS_* default and no pins; perturbations are root nodes reaching them (at most MAX_ROOTS, by stable
#       id), readouts are those nodes and up to MAX_TERMINALS terminal nodes they
#       reach. Every arm is then scored on the same pairs.
#   julia behaviour.jl BUILD run SELECTION.tsv OUT.tsv
#       solve every (pathway, root, dose) of the selection under the current DS_*
#       environment; one row per readout, plus the solve's convergence.
#   julia behaviour.jl BUILD relabel SELECTION.tsv OUT.tsv
#       per pathway, the first root at 2x solved on the network and on a copy with
#       every uuid renamed; the largest difference over all nodes.
using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal

const MAX_ROOTS = 8
const MAX_TERMINALS = 10
const DOSES = (0.0, 0.5, 2.0, 10.0)
params() = DS.SteadyStateParams(1.0, 0.1, 5000, 1e-6, "penalty")

build, mode = ARGS[1], ARGS[2]
load(pw) = (d = joinpath(build, pw);
            DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv")))
sid(net, u) = something(net.nodes[u].reactome_id, u)

function reach(adj, from)
    seen = Set(from); st = collect(from)
    while !isempty(st)
        v = pop!(st)
        for w in get(adj, v, String[]); w in seen || (push!(seen, w); push!(st, w)); end
    end
    seen
end

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

if mode == "select"
    out = DataFrame(pathway = String[], role = String[], uuid = String[], stid = String[])
    for pw in sort(filter(p -> startswith(p, "R-HSA-"), readdir(build)))
        isfile(joinpath(build, pw, "logic_network.csv")) || continue
        net = redirect_stdout(devnull) do; load(pw); end
        sn = DS.solved_carrier_network(first(DS.apply_material_loops(net)))
        held = DS.held_uuids(net)
        types = Dict{Tuple{String, String}, Set{String}}()
        sadj = Dict{String, Vector{String}}()
        for e in sn.edges
            push!(get!(types, (e.parent_uuid, e.child_uuid), Set{String}()), e.edge_type)
            (e.edge_type == DS.RECYCLE_EDGE || e.child_uuid in held) && continue
            push!(get!(sadj, e.parent_uuid, String[]), e.child_uuid)
        end
        selfmod = Set(p for (p, t) in types if "input" in t && "catalyst" in t)
        touched = Set{String}()
        for c in sccs(sort(collect(keys(sn.nodes))), sadj)
            cs = Set(c)
            any(p -> p[1] in cs && p[2] in cs, selfmod) && union!(touched, c)
        end
        isempty(touched) && continue
        adj = Dict{String, Vector{String}}(); radj = Dict{String, Vector{String}}()
        for e in net.edges
            e.edge_type == DS.RECYCLE_EDGE && continue
            push!(get!(adj, e.parent_uuid, String[]), e.child_uuid)
            push!(get!(radj, e.child_uuid, String[]), e.parent_uuid)
        end
        kind(u) = net.node_leaves === nothing || !haskey(net.node_leaves, u) ? "" : net.node_leaves[u][1]
        up = reach(radj, touched)
        roots = [u for u in up if !haskey(radj, u) && !(u in held) && kind(u) != "reaction"]
        sort!(roots; by = u -> (sid(net, u), u))
        down = reach(adj, touched)
        terms = [u for u in down if !haskey(adj, u) && !(u in held) && kind(u) != "reaction" && !(u in touched)]
        sort!(terms; by = u -> (sid(net, u), u))
        for u in roots[1:min(end, MAX_ROOTS)]; push!(out, (pw, "root", u, sid(net, u))); end
        for u in sort!(collect(touched); by = u -> (sid(net, u), u)); push!(out, (pw, "touched", u, sid(net, u))); end
        for u in terms[1:min(end, MAX_TERMINALS)]; push!(out, (pw, "terminal", u, sid(net, u))); end
        println(pw, ": touched ", length(touched), ", roots ", min(length(roots), MAX_ROOTS), " of ", length(roots),
                ", terminals ", min(length(terms), MAX_TERMINALS), " of ", length(terms))
    end
    CSV.write(ARGS[3], out; delim = '\t')
elseif mode in ("run", "relabel")
    sel = CSV.read(ARGS[3], DataFrame; delim = '\t', types = String)
    rows = DataFrame(pathway = String[], root = String[], dose = Float64[], readout = String[], role = String[],
                     fold = Float64[], converged = Bool[], residual = Float64[])
    for g in groupby(sel, :pathway)
        pw = g.pathway[1]
        net = redirect_stdout(devnull) do; load(pw); end
        roots = g.uuid[g.role .== "root"]
        reads = [(r.uuid, r.role) for r in eachrow(g) if r.role != "root"]
        isempty(roots) && continue
        if mode == "run"
            for (root, dose) in vcat([("", 1.0)], [(r, d) for r in roots for d in DOSES])
                obs = isempty(root) ? Dict{String, Tuple{Float64, Float64}}() : Dict(root => (dose, 1.0))
                res = DS.solve_steady_state(net, obs, params())
                for (u, role) in reads
                    push!(rows, (pw, root, dose, u, role, res.node_activities[u] / net.nodes[u].baseline,
                                 res.converged, res.final_residual))
                end
            end
        else
            root = roots[1]
            a = DS.solve_steady_state(net, Dict(root => (2.0, 1.0)), params())
            ren = Dict(u => "x" * reverse(u) for u in keys(net.nodes))
            rl(u) = get(ren, u, u)
            rn = DS.ReactionNetwork(
                Dict(rl(u) => DS.NetworkNode(rl(u), n.reactome_id, n.entity_type, n.original_set_id, n.display_name, n.baseline)
                     for (u, n) in net.nodes),
                [DS.LogicNetworkEdge(rl(e.parent_uuid), rl(e.child_uuid), e.is_and, e.is_positive, e.stoichiometry, e.edge_type)
                 for e in reverse(net.edges)],
                net.set_mappings, net.cofactor_stids, net.containment, net.drug_stids,
                net.pools === nothing ? nothing :
                    [DS.CyclePool(p.id, rl.(p.forms), rl(p.base),
                                  [DS.PoolPath(rl(q.from), rl(q.to), [(rl(a), rl(b), rl.(x)) for (a, b, x) in q.steps], q.enzyme)
                                   for q in p.paths],
                                  rl.(p.intermediates), [(rl(c), rl(x)) for (c, x) in p.carriers], rl.(p.bases)) for p in net.pools],
                net.pathogen_stids,
                net.carrier_loops === nothing ? nothing : [(rl(f), rl(s), c) for (f, s, c) in net.carrier_loops],
                net.node_leaves === nothing ? nothing : Dict(rl(u) => v for (u, v) in net.node_leaves),
                net.rna_stids, net.leaf_refs)
            b = DS.solve_steady_state(rn, Dict(rl(root) => (2.0, 1.0)), params())
            worst = maximum(abs(a.node_activities[u] - b.node_activities[rl(u)]) / net.nodes[u].baseline
                            for u in keys(net.nodes))
            push!(rows, (pw, root, 2.0, "max_abs_fold_difference", "relabel", worst, a.converged && b.converged,
                         max(a.final_residual, b.final_residual)))
        end
        println(pw, " done")
    end
    CSV.write(ARGS[4], rows; delim = '\t')
else
    error("mode must be select, run or relabel")
end
