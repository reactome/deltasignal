#!/usr/bin/env julia
# specs/142 parity gate, engine side: replay every fold of bench_rails.json
# through the engine's SERVED rail test (readout_railed with the env-resolved
# default stops — the exact function readouts_handler calls), and write the
# engine's flags beside the bench's for compare_142.py.
#
# Run in the dev container (the 141 pattern):
#   docker compose -f docker-compose.dev.yml run --rm \
#     -v $PWD/bench/analysis/142:/gate julia-api \
#     julia --project=/app /app/bench/analysis/142/rails_142.jl \
#       /gate/bench_rails.json /gate/engine_rails.json

using JSON3

include(joinpath(@__DIR__, "..", "..", "..", "src", "DeltaSignal.jl"))
using .DeltaSignal

function main()
    bench_path, out_path = ARGS[1], ARGS[2]
    bench = JSON3.read(read(bench_path, String))
    stops = DeltaSignal.readout_rail_stops()
    cases = [Dict(
        "axis" => c[:axis],
        "pathway" => c[:pathway],
        "gene" => c[:gene],
        "direction" => c[:direction],
        "key_output" => c[:key_output],
        "fold" => Float64(c[:fold]),
        "railed" => DeltaSignal.readout_railed(Float64(c[:fold]); stops = stops),
    ) for c in bench[:cases]]
    open(out_path, "w") do io
        JSON3.write(io, Dict("rail_low" => stops[1], "rail_high" => stops[2],
                             "cases" => cases))
    end
    println("$(length(cases)) folds flagged (stops $(stops)) -> $(out_path)")
end

main()
