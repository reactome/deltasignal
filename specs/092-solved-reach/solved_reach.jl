# specs/092: is a `no_path` case routed on the network the SOLVER runs?
#
# The bench labels a wrong case `no_path` when no edge path in the generator's
# logic network (plus boundary edges) leads from a pinned uuid to a readout uuid
# (bench/benchmark_vs_mpbiopath.py build_adjacency). The solver runs another
# network: recycle edges stripped, carrier-loop edges dropped, root twins,
# recycled supply, component twins and sequestration edges added
# (`solved_carrier_network`, the network influence scores use since the
# 2026-10-22 review). For every wrong case labelled no_path, and for a check on
# every wrong routed case, this asks whether a directed path exists on each.
#
# Usage (in the dev container, catalog mounted at /app/pathway_catalog):
#   julia --project=/app solved_reach.jl CASES_DIR OUT.tsv
# CASES_DIR holds curator_cases.tsv and experimental_cases.tsv.
using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal

cases_dir, out_path = ARGS[1], ARGS[2]
pw_ids = Dict(r.pathway_name => r.id for r in eachrow(CSV.read("/app/bench/catalog_pathways.tsv", DataFrame;
                                                                comment = "#", delim = '\t')))

function reach(edges, from::Vector{String})
    adj = Dict{String, Vector{String}}()
    for (a, b) in edges
        push!(get!(adj, a, String[]), b)
    end
    seen = Set(from); stack = copy(from)
    while !isempty(stack)
        v = pop!(stack)
        for w in get(adj, v, String[])
            w in seen || (push!(seen, w); push!(stack, w))
        end
    end
    seen
end


# The solver's pool update moves a protein along each path's steps (a -> copy ->
# b), and the last step lands in the target state with no logic-network edge.
pool_edges(net) = net.pools === nothing ? Tuple{String, String}[] :
    unique(vcat([[(a, rx) for (a, b, rxs) in q.steps for rx in rxs] for p in net.pools for q in p.paths]...,
                [[(rx, b) for (a, b, rxs) in q.steps for rx in rxs] for p in net.pools for q in p.paths]...))

rows = DataFrame[]
for ax in ("curator", "experimental")
    df = CSV.read(joinpath(cases_dir, "$(ax)_cases.tsv"), DataFrame; delim = '\t', types = String)
    df = df[(df.valid .== "1") .& (df.predicted .!= df.expected), :]
    df.axis .= ax
    push!(rows, df)
end
cases = vcat(rows...)
out = DataFrame(axis = String[], pathway = String[], gene = String[], direction = String[],
                key_output = String[], category = String[], raw = Bool[], solved = Bool[],
                solved_only_via = String[])
for pw in sort(unique(cases.pathway))
    haskey(pw_ids, pw) || (println(stderr, "no id for ", pw); continue)
    d = joinpath("/app/pathway_catalog", pw_ids[pw])
    isdir(d) || continue
    net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
    raw_edges = [(e.parent_uuid, e.child_uuid) for e in net.edges]          # as the bench reads them
    sn = DS.solved_carrier_network(net)
    solved_edges = vcat([(e.parent_uuid, e.child_uuid) for e in sn.edges], pool_edges(net))
    raw_set = Set(raw_edges)
    added = Dict((e.parent_uuid, e.child_uuid) => e.edge_type for e in sn.edges if !((e.parent_uuid, e.child_uuid) in raw_set))
    sub = cases[cases.pathway .== pw, :]
    memo = Dict{String, Tuple{Set{String}, Set{String}}}()
    for r in eachrow(sub)
        ismissing(r.gene_uuids) && continue
        genes = String.(split(r.gene_uuids, '|')); outs = String.(split(coalesce(r.output_uuids, ""), '|'))
        rr, rs = get!(memo, r.gene_uuids) do
            (reach(raw_edges, genes), reach(solved_edges, genes))
        end
        a = any(in(rr), outs); b = any(in(rs), outs)
        via = ""
        if b && !a
            # which added edge types the solved route needs: reach again without each type
            types = sort(unique(values(added)))
            need = [t for t in types if !any(in(reach([e for e in solved_edges if get(added, e, "") != t], genes)), outs)]
            via = join(need, "|")
        end
        push!(out, (r.axis, pw, r.gene, r.direction, r.key_output, r.category, a, b, via))
    end
end
CSV.write(out_path, out; delim = '\t')
println("cases ", nrow(out))
for ax in ("curator", "experimental")
    o = out[out.axis .== ax, :]
    np = o[o.category .== "no_path", :]
    println(ax, ": no_path ", nrow(np), "; routed on the generator network anyway ", count(np.raw),
            "; routed on the solved network ", count(np.solved))
    rt = o[o.category .!= "no_path", :]
    println(ax, ": routed-category wrong ", nrow(rt), "; with no route on the solved network ", count(.!rt.solved))
end

# Part 2: wrong cases routed on the generator network but not on the solved one.
# Which removal cuts them: restore one class of removed edges at a time.
cut = DataFrame(axis = String[], pathway = String[], gene = String[], direction = String[],
                key_output = String[], expected = String[], predicted = String[], restored_by = String[])
for pw in sort(unique(cases.pathway))
    haskey(pw_ids, pw) || continue
    d = joinpath("/app/pathway_catalog", pw_ids[pw]); isdir(d) || continue
    sub = cases[cases.pathway .== pw, :]
    net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
    sn = DS.solved_carrier_network(net)
    solved_edges = vcat([(e.parent_uuid, e.child_uuid) for e in sn.edges], pool_edges(net))
    sset = Set(solved_edges)
    removed = Dict{String, Vector{Tuple{String, String}}}()
    loop_pairs = net.carrier_loops === nothing ? Set{Tuple{String, String}}() :
                 Set((src, free) for (free, src, _) in net.carrier_loops)
    for e in net.edges
        k = (e.parent_uuid, e.child_uuid)
        k in sset && continue
        cls = e.edge_type == "recycle" ? "recycle" : k in loop_pairs ? "carrier_loop" : "other:" * e.edge_type
        push!(get!(removed, cls, Tuple{String, String}[]), k)
    end
    for r in eachrow(sub)
        ismissing(r.gene_uuids) && continue
        genes = String.(split(r.gene_uuids, '|')); outs = String.(split(coalesce(r.output_uuids, ""), '|'))
        any(in(reach([(e.parent_uuid, e.child_uuid) for e in net.edges], genes)), outs) || continue
        any(in(reach(solved_edges, genes)), outs) && continue
        by = [c for c in sort(collect(keys(removed))) if any(in(reach(vcat(solved_edges, removed[c]), genes)), outs)]
        push!(cut, (r.axis, pw, r.gene, r.direction, r.key_output, r.expected, r.predicted,
                    isempty(by) ? "only together" : join(by, "|")))
    end
end
CSV.write(replace(out_path, ".tsv" => "_cut.tsv"), cut; delim = '\t')
for ax in ("curator", "experimental")
    c = cut[cut.axis .== ax, :]
    println(ax, " cut by the solver: ", nrow(c))
    for (k, n) in sort(collect(Dict(k => count(==(k), c.restored_by) for k in unique(c.restored_by))); by = x -> -x[2])
        println("   restored by ", k, ": ", n)
    end
    for (k, n) in sort(collect(Dict(k => count(==(k), c.pathway) for k in unique(c.pathway))); by = x -> -x[2])[1:min(end, 8)]
        println("   ", k, ": ", n)
    end
end
