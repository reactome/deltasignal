# specs/095 census: inhibitor slots the copy clause adds, per pathway, with the
# solver's own pair test (self_inhibitor_setup under DS_SELF_INHIBITOR_COPIES=on).
# Usage (dev container, catalog at /app/pathway_catalog): julia census.jl
using DeltaSignal, CSV, DataFrames
const DS = DeltaSignal
tot = 0; rows = Tuple{String, Int, Int}[]
for pw in sort(filter(p -> startswith(p, "R-HSA-"), readdir("/app/pathway_catalog")))
    d = joinpath("/app/pathway_catalog", pw)
    isfile(joinpath(d, "logic_network.csv")) || continue
    net = DS.parse_complete_network(joinpath(d, "logic_network.csv"), joinpath(d, "stid_to_uuid_mapping.csv"))
    base = withenv(() -> DS.self_inhibitor_setup(net), "DS_SELF_INHIBITOR_COPIES" => "off")
    withc = withenv(() -> DS.self_inhibitor_setup(net), "DS_SELF_INHIBITOR_COPIES" => "on")
    n_new = count(k -> !haskey(base[1], k), keys(withc[1]))
    withc[4] > 0 && push!(rows, (pw, withc[4], n_new))
    global tot += withc[4]
end
for r in sort(rows; by = r -> -r[2]); println(r[1], "\t", r[2], " slots changed\t", r[3], " new slots"); end
println("total slots changed: ", tot, " in ", length(rows), " pathways")
