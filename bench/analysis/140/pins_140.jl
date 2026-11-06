#!/usr/bin/env julia
#
# specs/140 FR-003, engine side: resolve the SAME queries dump_pins.py dumped,
# through the parsed-network path the API serves (parse_complete_network +
# build_pin_index + resolve_pins), and write them for compare_140.py to diff.
#
# Run in the dev container, with the catalog the dump used mounted:
#   docker compose -f docker-compose.dev.yml run --rm -T julia-api \
#     julia --project=/app /app/bench/analysis/140/pins_140.jl \
#     /app/bench/analysis/140/bench_pins.json /app/bench/analysis/140/engine_pins.json

using JSON3

include(joinpath(@__DIR__, "..", "..", "..", "src", "DeltaSignal.jl"))
using .DeltaSignal

function main()
    in_path = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "bench_pins.json")
    out_path = length(ARGS) >= 2 ? ARGS[2] : joinpath(@__DIR__, "engine_pins.json")
    catalog = get(ENV, "DS_PATHWAY_CATALOG", "/app/pathway_catalog")
    dumped = JSON3.read(read(in_path, String))
    root_ignore = DeltaSignal.pin_root_ignore()
    println("catalog: $catalog  root_ignore: $(isempty(root_ignore) ? "none" : root_ignore)")

    out = []
    for p in dumped
        dir = joinpath(catalog, String(p[:dir]))
        isdir(dir) || error("$(p[:dir]) is not in $catalog; the engine must resolve the SAME build the dump used")
        network = DeltaSignal.parse_complete_network(
            joinpath(dir, "logic_network.csv"),
            joinpath(dir, "stid_to_uuid_mapping.csv"),
            nothing)
        index = DeltaSignal.build_pin_index(joinpath(dir, "nodes.csv"))
        index === nothing && error("$(p[:dir]) has no nodes.csv; the parity run needs the same build as the dump")
        queries = []
        for q in p[:queries]
            stids = String[String(s) for s in q[:stable_ids]]
            r = DeltaSignal.resolve_pins(stids, index, network.edges; root_ignore = root_ignore)
            push!(queries, Dict("gene" => String(q[:gene]),
                                "pins" => sort(unique(r.pins)),
                                "basis" => r.basis,
                                "occurrences" => r.occurrences))
        end
        push!(out, Dict("dir" => String(p[:dir]), "pathway_name" => String(p[:pathway_name]),
                        "queries" => queries))
        println("  $(p[:pathway_name]): $(length(queries)) genes")
    end
    open(out_path, "w") do io
        JSON3.write(io, out)
    end
    println("wrote $out_path: $(length(out)) pathways")
end

main()
