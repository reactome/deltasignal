#!/usr/bin/env julia
# specs/146 FR-004, engine side: replay every validation lookup through the
# SERVED readout_confidence (the shipped table) and write the answers.
# Run in the dev container (the 142 pattern), mounting this directory.
using JSON3
include(joinpath(@__DIR__, "..", "..", "..", "src", "DeltaSignal.jl"))
using .DeltaSignal

table, reason = DeltaSignal.load_confidence_table()
table === nothing && error("confidence table unreadable: $reason")
MAGFOLD = Dict("baseline" => 1.0, "le1" => 2.0, "le3" => 4.0, "gt3" => 100.0)
cases = JSON3.read(read(ARGS[1], String))
out = [begin
    # a representative fold inside the case's bin reproduces the bin
    conf, stratum = DeltaSignal.readout_confidence(table, String(c[:call]), Bool(c[:railed]),
                                                   String(c[:basis]), MAGFOLD[String(c[:mag])];
                                                   converged = Bool(c[:conv]))
    Dict("confidence" => conf === nothing ? nothing : round(conf, digits = 6),
         "stratum" => stratum)
end for c in cases]
open(ARGS[2], "w") do io; JSON3.write(io, out); end
println("$(length(out)) lookups replayed")
