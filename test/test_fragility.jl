#!/usr/bin/env julia
#
# specs/157: structural fragility — the table loader's honest refusals
# (absent, malformed, provenance mismatch, unverifiable serving record), the
# lookup's level walk on a fixture table whose cells are COMPUTED here from
# explicit fixture labels (never pinned), the minimum-cell refusal, the
# untouched rule, and the SHIPPED table's provenance pinned to the serving
# record. The gate measurement lives in the spec's record
# (bench/analysis/157, specs/157 research.md), not here.

using Test
using JSON3

include(joinpath(@__DIR__, "..", "src", "DeltaSignal.jl"))
using .DeltaSignal

const SERVED_BUILD = "20261113-1311_baf608b_r98fix"

# Fixture labels the test can see: 10 pairs in the down|le3|True class, 3 of
# them fragile — the cell's rate is computed from these, never asserted twice.
fixture_labels = [true, true, true, false, false, false, false, false, false, false]
RATE = count(fixture_labels) / length(fixture_labels)

function fixture_table(; min_cell = 5, catalog_build = SERVED_BUILD, drop = String[])
    cell = Dict("fragile_rate" => RATE, "n" => length(fixture_labels))
    small = Dict("fragile_rate" => 1.0, "n" => 2)
    t = Dict{String, Any}(
        "version" => "fixture", "spec" => "specs/157-structural-fragility",
        "catalog_build" => catalog_build, "min_cell" => min_cell,
        "levels" => ["full", "no_converged", "call_only"],
        "table" => Dict(
            # the full cell is TOO SMALL, so a converged-carrying lookup must
            # fall through to no_converged — the walk is what is tested
            "full" => Dict("down|le3|True" => small),
            "no_converged" => Dict("down|le3" => cell),
            "call_only" => Dict("down" => cell)))
    for k in drop
        delete!(t, k)
    end
    t
end

roundtrip(t) = JSON3.read(JSON3.write(t))

function with_table(f, t)
    path = tempname() * ".json"
    write(path, JSON3.write(t))
    try
        f(path)
    finally
        rm(path; force = true)
    end
end

@testset "structural fragility (specs/157)" begin

    @testset "loader refusals are honest and specific" begin
        t, reason = DeltaSignal.load_fragility_table(tempname() * ".json";
                                                     served_build = SERVED_BUILD)
        @test t === nothing && reason == "missing"

        with_table(Dict("not" => "a table")) do path
            t, reason = DeltaSignal.load_fragility_table(path; served_build = SERVED_BUILD)
            @test t === nothing && reason == "malformed"
        end

        with_table(fixture_table(catalog_build = "someone-elses-build")) do path
            t, reason = DeltaSignal.load_fragility_table(path; served_build = SERVED_BUILD)
            @test t === nothing
            @test occursin("provenance mismatch", reason)
            @test occursin("someone-elses-build", reason)
        end

        with_table(fixture_table()) do path
            t, reason = DeltaSignal.load_fragility_table(path; served_build = nothing)
            @test t === nothing
            @test occursin("unrecorded", reason)
            t2, reason2 = DeltaSignal.load_fragility_table(path; served_build = "unrecorded")
            @test t2 === nothing
        end

        with_table(fixture_table()) do path
            t, reason = DeltaSignal.load_fragility_table(path; served_build = SERVED_BUILD)
            @test t !== nothing && reason == ""
        end
    end

    @testset "the level walk, computed from the fixture labels" begin
        t = roundtrip(fixture_table())
        # full cell (n=2) is below min_cell=5: a converged lookup falls
        # through to no_converged and serves the computed rate
        fr = DeltaSignal.readout_fragility(t, "down", 0.4; converged = true)
        @test fr !== nothing
        @test fr[1] ≈ RATE
        @test fr[2] == length(fixture_labels)
        # without converged the walk starts at no_converged: same cell
        @test DeltaSignal.readout_fragility(t, "down", 0.4) == fr
        # the endpoint's "normal" vocabulary maps to the record's "not":
        # no cell for it in this fixture, so honest absence
        @test DeltaSignal.readout_fragility(t, "normal", 0.4) === nothing
        # a magnitude outside the cell's bin falls to call_only (same rate here)
        fr2 = DeltaSignal.readout_fragility(t, "down", 0.05; converged = true)
        @test fr2 !== nothing && fr2[1] ≈ RATE
    end

    @testset "minimum cell and the untouched rule" begin
        t = roundtrip(fixture_table(min_cell = 50))   # every cell too small
        @test DeltaSignal.readout_fragility(t, "down", 0.4; converged = true) === nothing
        t2 = roundtrip(fixture_table())
        # fold exactly 1.0: untouched, nothing served whatever the table holds
        @test DeltaSignal.readout_fragility(t2, "down", 1.0; converged = true) === nothing
    end

    @testset "the shipped table is the serving record's" begin
        shipped, reason = DeltaSignal.load_fragility_table(; served_build = SERVED_BUILD)
        @test shipped !== nothing && reason == ""
        @test string(shipped[:spec]) == "specs/157-structural-fragility"
        @test string(shipped[:catalog_build]) == SERVED_BUILD
        min_cell = Int(shipped[:min_cell])
        @test min_cell == 30
        n_cells = 0
        for level in keys(shipped[:table]), (_, cell) in shipped[:table][level]
            @test Int(cell[:n]) >= min_cell
            @test 0.0 <= Float64(cell[:fragile_rate]) <= 1.0
            n_cells += 1
        end
        @test n_cells > 0
        # the measured headline classes exist: every call serves at call_only
        for call in ("down", "not", "up")
            fr = DeltaSignal.readout_fragility(shipped, call, 0.4)
            call == "down" && @test fr !== nothing
        end
        # and a mismatch against a different build refuses the SHIPPED table
        t, r = DeltaSignal.load_fragility_table(; served_build = "not-this-build")
        @test t === nothing && occursin("provenance mismatch", r)
    end
end

println("test_fragility.jl passed")
