using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal
d = "/app/pathway_catalog/R-HSA-1257604"
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
nd = CSV.read(joinpath(d, "nodes.csv"), DataFrame); ent = Dict(r.uuid => string(r.diagram_entity_id) for r in eachrow(nd))
pins, outs = split.(split(strip(read("/scratch/tp53pip3.tsv", String)), '\t'), '|')
println("TP53 pins: ", [ent[p] for p in pins])
byp(p) = first(u for u in keys(net.nodes) if startswith(u, p))
watch = [("PTEN mRNA root", byp("08cc6d72")), ("PTEN mRNA transcribed", byp("fbc1de1c")), ("miR complex", byp("a04983dd")),
         ("translation", byp("9cc87071")), ("PTEN copy (acts on PIP3)", byp("57a5dff3")), ("PIP3", byp("d7b2a6f6"))]
for v in (1.0, 80.0, 0.0)
    b = DS.solve_steady_state(net, Dict(String(p) => (v, 1.0) for p in pins))
    f(u) = round(b.node_activities[u] / 0.01, digits = 3)
    println("TP53=$v: ", [(k, f(u)) for (k, u) in watch], " readout ", [f(u) for u in outs])
end
