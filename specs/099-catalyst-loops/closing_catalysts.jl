# specs/099: which catalyst edges close the unhandled catalyst-closed positive loops
# (specs/097), and how the catalyst relates to the reaction it catalyses.
#
# For each unhandled positive component whose cycles all need a catalyst or
# regulator edge, every catalyst or positive-regulator edge E -> R inside it (each lies on a cycle,
# since both ends are in one strongly connected component) is classified:
#   self      - E is also an input of R (the complex is curated as its own catalyst)
#   regen     - E's stable id is among R's outputs (released unchanged)
#   contains  - an output of R contains E's entity, or E contains an input's (containment)
#   feedback  - none of these: R's products lead to a different entity that activates E
# Usage (dev container): julia closing_catalysts.jl BUILD_DIR OUT.tsv
using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal
build, out_path = ARGS[1], ARGS[2]
include(joinpath(@__DIR__, "..", "097-loop-audit", "sccs.jl"))
rows = DataFrame(pathway = String[], comp_size = Int[], role = String[], catalyst = String[], reaction = String[], kind = String[])
for pw in sort(filter(p -> startswith(p, "R-HSA-"), readdir(build)))
    d = joinpath(build, pw); isfile(joinpath(d, "logic_network.csv")) || continue
    net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
    sn = DS.solved_carrier_network(net)
    sid(u) = (n = get(sn.nodes, u, nothing); n === nothing ? "" : something(n.reactome_id, ""))
    pool_of = Dict{String, Int}()
    if net.pools !== nothing
        for (k, p) in enumerate(net.pools)
            for u in vcat(p.forms, p.intermediates); pool_of[u] = k; end
            for q in p.paths, (a, b, rxs) in q.steps, rx in rxs; pool_of[rx] = k; end
        end
    end
    self_inh = first(DS.self_inhibitor_setup(sn)); self_act = first(DS.self_activator_setup(sn))
    adj = Dict{String, Vector{String}}(); mat = Dict{String, Vector{String}}()
    neg = Set{Tuple{String, String}}(); et = Dict{Tuple{String, String}, String}()
    ins = Dict{String, Vector{String}}(); outs = Dict{String, Vector{String}}()
    for e in sn.edges
        a, b = e.parent_uuid, e.child_uuid
        e.edge_type == "input" && push!(get!(ins, b, String[]), a)
        e.edge_type == "output" && push!(get!(outs, a, String[]), b)
        (haskey(pool_of, a) && get(pool_of, b, 0) == pool_of[a]) && continue
        ((!e.is_positive && haskey(self_inh, (a, b))) || (e.is_positive && haskey(self_act, (a, b)))) && continue
        push!(get!(adj, a, String[]), b); et[(a, b)] = e.edge_type
        e.is_positive || push!(neg, (a, b))
        e.edge_type in ("catalyst", "regulator") || push!(get!(mat, a, String[]), b)
    end
    cont = sn.containment
    for c in sccs(sort(collect(keys(sn.nodes))), adj)
        cs = Set(c)
        inner = [(a, b) for a in c for b in get(adj, a, String[]) if b in cs]
        any(in(neg), inner) && continue                    # positive components only
        madj = Dict(a => [b for b in get(mat, a, String[]) if b in cs] for a in c)
        isempty(sccs(sort(c), madj)) || continue           # catalyst-closed only
        for (E, R) in inner
            et[(E, R)] in ("catalyst", "regulator") || continue
            es = sid(E); is = [sid(i) for i in get(ins, R, String[])]; os = [sid(o) for o in get(outs, R, String[])]
            kind = E in get(ins, R, String[]) || es in is ? "self" :
                   es in os ? "regen" :
                   (any(o -> es in get(cont, o, Set{String}()), os) || any(i -> i in get(cont, es, Set{String}()), is)) ? "contains" :
                   "feedback"
            push!(rows, (pw, length(c), et[(E, R)], es, sid(R), kind))
        end
    end
end
CSV.write(out_path, rows; delim = '\t')
println("catalyst and positive-regulator edges inside catalyst-closed components: ", nrow(rows), " in ", length(unique(rows.pathway)), " pathways")
for r in ("catalyst", "regulator"), k in ("self", "regen", "contains", "feedback")
    x = rows[(rows.kind .== k) .& (rows.role .== r), :]
    nrow(x) == 0 && continue
    println("  ", rpad(r, 10), rpad(k, 9), nrow(x), " edges, ", length(unique(x.pathway)), " pathways, ",
            length(unique(zip(x.catalyst, x.reaction))), " distinct")
end
println("pathways with any: ", length(unique(rows.pathway)))
