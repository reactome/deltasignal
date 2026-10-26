# specs/107 amendment 1: rule 3b's Class I MHC rows (UBA7, BAP31) added one at a time; the case uuids come from
# cblb.txt (gene, readout, pinned uuids, readout uuids, from mc107's curator_cases.tsv). Mount the refs107m build at /build.
using DeltaSignal, Statistics; const DS = DeltaSignal
d = "/build/R-HSA-983169"
net = redirect_stdout(devnull) do; DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv")); end
rows, _ = DS.multistep_carriers(net, "compartment")
with(rs) = DS.ReactionNetwork(net.nodes, net.edges, net.set_mappings, net.cofactor_stids, net.containment, net.drug_stids,
                              net.pools, net.pathogen_stids, vcat(something(net.carrier_loops, DS.CarrierLoopEdge[]), rs),
                              net.node_leaves, net.rna_stids, net.leaf_refs)
uba7 = [r for r in rows if net.nodes[r[1]].reactome_id == "R-HSA-936550"]
bap = [r for r in rows if net.nodes[r[1]].reactome_id == "R-HSA-983323"]
ENV["DS_CARRIER_MULTISTEP"] = "off"
for line in eachline("/scratch/cblb.txt")
    g, ko, pins, outs = split(line)
    obs = Dict(String(u) => (0.0, 1.0) for u in split(pins, "|"))
    o = split(outs, "|")
    for (lab, rs) in (("none", DS.CarrierLoopEdge[]), ("UBA7", uba7), ("BAP31", bap), ("both", rows))
        r = DS.solve_steady_state(with(rs), obs, DS.SteadyStateParams(1.0, 0.1, 5000, 1e-6, "penalty"))
        println(g, " ", ko, " ", rpad(lab, 6), " mean fold ", round(mean(r.node_activities[u] for u in o) / 0.01, digits = 4), " conv ", r.converged)
    end
end
