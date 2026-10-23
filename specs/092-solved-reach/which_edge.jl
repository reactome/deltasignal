# specs/092: for one wrong case cut by the solver, which single dropped carrier-loop
# edge restores its route. Usage: julia which_edge.jl PATHWAY_ID GENE_UUIDS OUTPUT_UUIDS
using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal
d = joinpath("/app/pathway_catalog", ARGS[1])
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
nd = CSV.read(joinpath(d, "nodes.csv"), DataFrame); ent = Dict(r.uuid => string(r.diagram_entity_id) for r in eachrow(nd))
genes = String.(split(ARGS[2], '|')); outs = String.(split(ARGS[3], '|'))
pe = net.pools === nothing ? Tuple{String,String}[] :
    unique(vcat([[(a, rx) for (a, b, rxs) in q.steps for rx in rxs] for p in net.pools for q in p.paths]...,
                [[(rx, b) for (a, b, rxs) in q.steps for rx in rxs] for p in net.pools for q in p.paths]...))
sn = DS.solved_carrier_network(net)
E = vcat([(e.parent_uuid, e.child_uuid) for e in sn.edges], pe)
function reach(edges)
    adj = Dict{String, Vector{String}}(); for (a, b) in edges; push!(get!(adj, a, String[]), b); end
    seen = Set(genes); st = copy(genes)
    while !isempty(st); v = pop!(st); for w in get(adj, v, String[]); w in seen || (push!(seen, w); push!(st, w)); end; end
    seen
end
println("routed on solved: ", any(in(reach(E)), outs))
sset = Set(E)
drop = [(src, free, closed) for (free, src, closed) in net.carrier_loops if !((src, free) in sset)]
for (tag, sel) in (("closed", x -> x[3]), ("open", x -> !x[3]))
    add = [(a, b) for (a, b, c) in drop if sel((a, b, c))]
    println("restoring all ", tag, "-loop drops (", length(add), "): ", any(in(reach(vcat(E, add))), outs))
end
# by free form: which free form's dropped in-edges restore it
for free in unique(b for (a, b, c) in drop)
    add = [(a, b) for (a, b, c) in drop if b == free]
    any(in(reach(vcat(E, add))), outs) &&
        println("restores: all drops into ", ent[free], " (", length(add), " edges, ",
                any(c for (a, b, c) in drop if b == free) ? "closed" : "open", ")")
end
