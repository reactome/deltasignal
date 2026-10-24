using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal
d = "/app/pathway_catalog/R-HSA-453279"
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
nd = CSV.read(joinpath(d, "nodes.csv"), DataFrame); ent = Dict(r.uuid => string(r.diagram_entity_id) for r in eachrow(nd))
pins, outs = split.(split(strip(read("/scratch/ccne1.tsv", String)), '\t'), '|')
sn = DS.solved_carrier_network(net)
UB = Set(["R-HSA-68524", "R-HSA-113595", "R-HSA-9660007"])
fw = Dict{String, Vector{String}}(); bw = Dict{String, Vector{String}}()
for e in sn.edges; push!(get!(fw, e.parent_uuid, String[]), e.child_uuid); push!(get!(bw, e.child_uuid, String[]), e.parent_uuid); end
function reach(adj, from, block)
    seen = Set(from); st = collect(from)
    while !isempty(st); v = pop!(st); for w in get(adj, v, String[]); (w in seen || w in block) && continue; push!(seen, w); push!(st, w); end; end
    seen
end
ubn = Set(u for u in keys(sn.nodes) if ent[u] in UB)
f = reach(fw, String.(pins), Set{String}()); b = reach(bw, String.(outs), Set{String}())
println("route with free Ub: ", !isempty(intersect(f, b)))
f2 = reach(fw, String.(pins), ubn)
println("route without passing through free Ub: ", any(in(f2), outs))
on = intersect(f, b)
println("free Ub nodes on the route: ", [(ent[u], u[1:8]) for u in on if u in ubn])
# reactions producing / consuming those Ub nodes on the route
for u in on
    u in ubn || continue
    prod = [ent[e.parent_uuid] for e in sn.edges if e.child_uuid == u && e.parent_uuid in on]
    cons = [ent[e.child_uuid] for e in sn.edges if e.parent_uuid == u && e.child_uuid in on]
    println("  ", u[1:8], " made by ", unique(prod)[1:min(end,4)], " used by ", unique(cons)[1:min(end,4)])
end
