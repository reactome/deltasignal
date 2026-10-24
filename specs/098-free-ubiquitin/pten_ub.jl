using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal
d = "/app/pathway_catalog/R-HSA-1257604"
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
nd = CSV.read(joinpath(d, "nodes.csv"), DataFrame); ent = Dict(r.uuid => string(r.diagram_entity_id) for r in eachrow(nd))
pins, outs = split.(split(strip(read("/scratch/pten.tsv", String)), '\t'), '|')
byp(p) = first(u for u in keys(net.nodes) if startswith(u, p))
ub = sort([u for u in keys(net.nodes) if ent[u] in ("R-HSA-113595", "R-HSA-68524", "R-HSA-9660007")])
watch = [("translation", byp("9cc87071")), ("PTEN copy", byp("57a5dff3")), ("p-MKRN1", byp("45520a0c")), ("PIP3", byp("d7b2a6f6"))]
for m in ("propagate", "inert"), v in (1.0, 80.0)
    b = withenv("DS_UBIQUITIN_MODE" => m) do; DS.solve_steady_state(net, Dict(String(p) => (v, 1.0) for p in pins)); end
    f(u) = round(b.node_activities[u] / 0.01, digits = 3)
    println("ub=$m PTEN=$v conv=$(b.converged): ", [(k, f(u)) for (k, u) in watch], " readout ", [f(u) for u in outs],
            " free Ub ", [(ent[u], f(u)) for u in ub])
end
