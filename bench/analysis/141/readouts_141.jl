#!/usr/bin/env julia
#
# specs/141 FR-003/FR-004, engine side: resolve the SAME readouts
# dump_readouts.py dumped, through the path the API serves
# (parse_complete_network + build_pin_index + build_readout_index +
# resolve_readout / resolve_abundance_occurrences + abundance_copies), and
# write them for compare_141.py to diff. Set members come from the dump —
# the bench's own Neo4j answer — because the engine has no Neo4j by design.
#
# Run in the dev container, with the catalog the dump used mounted:
#   docker compose -f docker-compose.dev.yml run --rm -T julia-api \
#     julia --project=/app /app/bench/analysis/141/readouts_141.jl \
#     /app/bench/analysis/141/bench_readouts.json /tmp/engine_readouts.json

using JSON3

include(joinpath(@__DIR__, "..", "..", "..", "src", "DeltaSignal.jl"))
using .DeltaSignal

function main()
    in_path = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "bench_readouts.json")
    out_path = length(ARGS) >= 2 ? ARGS[2] : joinpath(@__DIR__, "engine_readouts.json")
    catalog = get(ENV, "DS_PATHWAY_CATALOG", "/app/pathway_catalog")
    dumped = JSON3.read(read(in_path, String))
    println(stderr, "catalog: $catalog")

    curator = []
    for p in dumped[:curator]
        dir = joinpath(catalog, String(p[:dir]))
        isdir(dir) || error("$(p[:dir]) is not in $catalog; the engine must resolve the SAME build the dump used")
        network = DeltaSignal.parse_complete_network(
            joinpath(dir, "logic_network.csv"),
            joinpath(dir, "stid_to_uuid_mapping.csv"),
            nothing)
        pin_index = DeltaSignal.build_pin_index(joinpath(dir, "nodes.csv"))
        pin_index === nothing && error("$(p[:dir]) has no nodes.csv; the parity run needs the same build as the dump")
        ridx = DeltaSignal.build_readout_index(dir)
        readouts = []
        for q in p[:readouts]
            stids = String[String(s) for s in q[:stable_ids]]
            members = String[String(m) for m in q[:members]]
            r = DeltaSignal.resolve_readout(stids, members, pin_index, ridx, network.edges)
            push!(readouts, Dict("ko" => String(q[:ko]),
                                 "basis" => r.basis,
                                 "groups" => r.groups))
        end
        push!(curator, Dict("dir" => String(p[:dir]), "readouts" => readouts))
        println(stderr, "  curator $(p[:dir]): $(length(readouts)) readouts")
    end

    korkut = []
    for p in dumped[:korkut]
        dir = joinpath(catalog, String(p[:dir]))
        isdir(dir) || error("$(p[:dir]) is not in $catalog")
        network = DeltaSignal.parse_complete_network(
            joinpath(dir, "logic_network.csv"),
            joinpath(dir, "stid_to_uuid_mapping.csv"),
            nothing)
        pin_index = DeltaSignal.build_pin_index(joinpath(dir, "nodes.csv"))
        pin_index === nothing && error("$(p[:dir]) has no nodes.csv")
        readouts = []
        for q in p[:readouts]
            stids = String[String(s) for s in q[:stable_ids]]
            occ = DeltaSignal.resolve_abundance_occurrences(stids, pin_index)
            kept = DeltaSignal.abundance_copies(occ, network.edges)
            push!(readouts, Dict("readout" => String(q[:readout]),
                                 "occurrences" => occ,
                                 "kept" => kept))
        end
        push!(korkut, Dict("dir" => String(p[:dir]), "readouts" => readouts))
        println(stderr, "  korkut $(p[:dir]): $(length(readouts)) total readouts")
    end

    open(out_path, "w") do io
        JSON3.write(io, Dict("curator" => curator, "korkut" => korkut))
    end
    println(stderr, "wrote $out_path: $(length(curator)) curator pathways, $(length(korkut)) korkut pathways")
end

main()
