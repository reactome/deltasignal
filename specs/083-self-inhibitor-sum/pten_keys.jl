using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal
d = "/app/pathway_catalog/R-HSA-1257604"
net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
pins = split(split(strip(read("/scratch/pten.tsv", String)), '\t')[1], '|')
nd = CSV.read(joinpath(d, "nodes.csv"), DataFrame); ent = Dict(r.uuid => string(r.diagram_entity_id) for r in eachrow(nd))
keys_ = ["b023efea" => "transcription", "94050adf" => "PTEN mRNA (made)", "bf7c5250" => "miR-26A:PTEN mRNA",
         "26483e7b" => "translation", "b728e14d" => "PTEN (made)", "63f70ea7" => "p-MKRN1", "82e92c08" => "PIP3", "dcdc49eb" => "p-AKT1"]
for v in (1.01, 2.0, 80.0)
    r = DS.solve_steady_state(net, Dict(String(p) => (v, 1.0) for p in pins))
    println("PTEN ", v, ": ", join([string(n, " ", round(r.node_activities[first(u for u in keys(net.nodes) if startswith(u, k))] / 0.01, digits = 3)) for (k, n) in keys_], " | "))
end
