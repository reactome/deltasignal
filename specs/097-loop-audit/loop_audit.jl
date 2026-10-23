# specs/097: loop audit. Which feedback loops does each network contain, and
# which of them does a solver rule handle?
#
# For each pathway in a build:
#   1. the network the solver runs (`solved_carrier_network`: recycle edges
#      stripped, carrier-loop edges dropped, root and component twins and
#      recycled supply added);
#   2. its non-trivial strongly connected components (cyclic components);
#   3. the edges a rule handles are removed:
#        pool  - an edge between two nodes of one interconversion, recycle or
#                release pool (states, intermediates, step copies; specs/039,
#                081, 086), whose cycle the pool balance solves;
#        self  - an inhibitor or activator slot built from its own reaction's
#                input (specs/022, 040, 083, 095), counted once by the rule;
#   4. the cyclic components left are UNHANDLED loops, labelled by shape:
#        positive      - every edge positive: gain 1 at baseline, the knife-edge
#                        (specs/014, RAF, specs/054);
#        self_depleting- holds a depletion edge from a complex onto a node of an
#                        entity it contains (a Goodwin ring, specs/066: the
#                        non-converged fallback finds its fixed point);
#        negative      - other negative feedback (an inhibitor or depletion edge).
#
# With a cases directory (a benchmark build), each case is assigned the unhandled
# loop types its route passes through (a node of an unhandled component that the
# perturbed nodes reach and that reaches a readout), or `no_route` (nothing reaches
# a readout on the solved network), or `no_unhandled_loop`.
#
# Usage (dev container): julia loop_audit.jl BUILD_DIR OUT_PREFIX [CASES_DIR]
using DeltaSignal, CSV, DataFrames, Statistics
const DS = DeltaSignal

build, prefix = ARGS[1], ARGS[2]
cases_dir = length(ARGS) >= 3 ? ARGS[3] : ""

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

pw_rows = DataFrame(pathway = String[], nodes = Int[], cyclic_nodes = Int[], cyclic_components = Int[],
                    handled_pool_edges = Int[], handled_self_edges = Int[],
                    unhandled_components = Int[], unhandled_nodes = Int[],
                    positive = Int[], self_depleting = Int[], negative = Int[])
loop_rows = DataFrame(pathway = String[], size = Int[], label = String[], neg_edges = Int[],
                      edge_types = String[], example = String[])
case_types = Dict{Tuple{String, String, String, String}, String}()   # (pathway dir, gene, dir, readout) -> types

cases = nothing
if !isempty(cases_dir)
    local cs = DataFrame[]
    for ax in ("curator", "experimental")
        f = joinpath(cases_dir, "$(ax)_cases.tsv"); isfile(f) || continue
        df = CSV.read(f, DataFrame; delim = '\t', types = String); df.axis .= ax
        push!(cs, df[df.valid .== "1", :])
    end
    cases = vcat(cs...)
    pw_ids = Dict(r.pathway_name => r.id for r in eachrow(CSV.read("/app/bench/catalog_pathways.tsv", DataFrame;
                                                                    comment = "#", delim = '\t')))
    cases.pid = [get(pw_ids, p, "") for p in cases.pathway]
end

for pw in sort(filter(p -> startswith(p, "R-HSA-"), readdir(build)))
    d = joinpath(build, pw)
    isfile(joinpath(d, "logic_network.csv")) || continue
    net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
    sn = DS.solved_carrier_network(net)
    sid(u) = (n = get(sn.nodes, u, nothing); n === nothing ? "" : something(n.reactome_id, ""))
    # pool node sets
    pool_of = Dict{String, Int}()
    if net.pools !== nothing
        for (k, p) in enumerate(net.pools)
            for u in vcat(p.forms, p.intermediates); pool_of[u] = k; end
            for q in p.paths, (a, b, rxs) in q.steps, rx in rxs; pool_of[rx] = k; end
        end
    end
    self_inh = first(DS.self_inhibitor_setup(sn))
    self_act = first(DS.self_activator_setup(sn))
    handled_pool = 0; handled_self = 0
    adj_all = Dict{String, Vector{String}}(); adj = Dict{String, Vector{String}}()
    neg = Set{Tuple{String, String}}(); etype = Dict{Tuple{String, String}, String}()
    for e in sn.edges
        a, b = e.parent_uuid, e.child_uuid
        push!(get!(adj_all, a, String[]), b)
        etype[(a, b)] = e.edge_type
        e.is_positive || push!(neg, (a, b))
        if haskey(pool_of, a) && get(pool_of, b, 0) == pool_of[a]
            handled_pool += 1; continue
        end
        if (!e.is_positive && haskey(self_inh, (a, b))) || (e.is_positive && haskey(self_act, (a, b)))
            handled_self += 1; continue
        end
        push!(get!(adj, a, String[]), b)
    end
    nodes = sort(collect(keys(sn.nodes)))
    comps_all = sccs(nodes, adj_all)
    comps = sccs(nodes, adj)
    cont = sn.containment
    labels = String[]
    for c in comps
        cs = Set(c)
        inner = [(a, b) for a in c for b in get(adj, a, String[]) if b in cs]
        nneg = count(in(neg), inner)
        selfdep = any(((a, b),) -> etype[(a, b)] == "depletion" &&
                      sid(b) in get(cont, sid(a), Set{String}()), inner)
        lab = nneg == 0 ? "positive" : selfdep ? "self_depleting" : "negative"
        push!(labels, lab)
        types = join(sort(unique([etype[x] for x in inner])), "|")
        push!(loop_rows, (pw, length(c), lab, nneg, types, sid(first(sort(c)))))
    end
    push!(pw_rows, (pw, length(nodes), sum(length, comps_all; init = 0), length(comps_all), handled_pool, handled_self,
                    length(comps), sum(length, comps; init = 0),
                    count(==("positive"), labels), count(==("self_depleting"), labels), count(==("negative"), labels)))
    if cases !== nothing
        sub = cases[cases.pid .== pw, :]
        comp_of = Dict(u => (k, labels[k]) for (k, c) in enumerate(comps) for u in c)
        radj = Dict{String, Vector{String}}()
        for (a, ws) in adj_all, b in ws; push!(get!(radj, b, String[]), a); end
        memo_f = Dict{String, Set{String}}(); memo_b = Dict{String, Set{String}}()
        for r in eachrow(sub)
            (ismissing(r.gene_uuids) || ismissing(r.output_uuids)) && continue
            fw = get!(() -> reach(adj_all, String.(split(r.gene_uuids, '|'))), memo_f, r.gene_uuids)
            bw = get!(() -> reach(radj, String.(split(r.output_uuids, '|'))), memo_b, r.output_uuids)
            on = intersect(fw, bw)
            ts = sort(unique([comp_of[u][2] for u in on if haskey(comp_of, u)]))
            # no route at all (the bench's no_path, on the solved network), else the
            # unhandled loop types on the route, else an acyclic or handled route
            case_types[(pw, r.gene, r.direction, r.key_output)] =
                isempty(on) ? "no_route" : isempty(ts) ? "no_unhandled_loop" : join(ts, "+")
        end
    end
end
CSV.write(prefix * "_pathways.tsv", pw_rows; delim = '\t')
CSV.write(prefix * "_loops.tsv", loop_rows; delim = '\t')
println("pathways ", nrow(pw_rows), "; cyclic components ", sum(pw_rows.cyclic_components),
        "; unhandled ", sum(pw_rows.unhandled_components), " (positive ", sum(pw_rows.positive),
        ", self_depleting ", sum(pw_rows.self_depleting), ", negative ", sum(pw_rows.negative), ")")
println("pathways with any unhandled loop: ", count(>(0), pw_rows.unhandled_components), " of ", nrow(pw_rows))
for lab in ("positive", "self_depleting", "negative")
    println("  pathways with a $lab unhandled loop: ", count(>(0), pw_rows[!, Symbol(lab)]), " of ", nrow(pw_rows))
end
if cases !== nothing
    cases.loop = [get(case_types, (r.pid, r.gene, r.direction, r.key_output), "unscored") for r in eachrow(cases)]
    out = combine(groupby(cases, [:axis, :loop]), nrow => :cases,
                  [:predicted, :expected] => ((p, e) -> round(mean(p .== e), digits = 4)) => :accuracy)
    sort!(out, [:axis, order(:cases, rev = true)])
    CSV.write(prefix * "_cases.tsv", out; delim = '\t')
    show(out; allrows = true); println()
end
