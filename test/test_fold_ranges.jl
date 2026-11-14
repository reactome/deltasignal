#!/usr/bin/env julia
#
# specs/155: measured fold intervals — the table loader's honest refusals
# (absent, malformed, provenance mismatch, unverifiable serving record), the
# lookup's stratum walk and arm-(b) band rule on a fixture table whose cells
# are COMPUTED here from fixture pairs (never pinned), and the SHIPPED
# table's provenance pinned to the serving record's constants. The
# out-of-sample gates and the end-to-end parity live in the spec's record
# (bench/analysis/155, specs/155 research.md), not here.

using Test
using JSON3

include(joinpath(@__DIR__, "..", "src", "DeltaSignal.jl"))
using .DeltaSignal

const SERVED_BUILD = "20261113-1311_baf608b_r98fix"
const SERVED_SOLVER = "12b5db3"

# A fixture table computed from explicit pairs: the down|False|le3|True cell
# holds measured folds 0.2..0.8 (middle 50% = [0.35, 0.65] by linear
# interpolation over the 7 sorted values), so the asserted interval is the
# quantile of data the test can see, not a number the table asserts twice.
fixture_measured = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
q(vals, p) = begin
    s = sort(vals); pos = p * (length(s) - 1)
    lo = Int(floor(pos)); hi = min(lo + 1, length(s) - 1)
    s[lo + 1] + (s[hi + 1] - s[lo + 1]) * (pos - lo)
end
cell_lo, cell_hi = q(fixture_measured, 0.25), q(fixture_measured, 0.75)

function fixture_table(; min_cell = 5, catalog_build = SERVED_BUILD,
                       arm_b = false, drop = String[])
    cell = Dict("n" => length(fixture_measured),
                "shares" => Dict("0.5" => Dict("lo" => cell_lo, "hi" => cell_hi)))
    if arm_b
        cell["bands"] = [Dict("fold_lo" => 0.1, "fold_hi" => 0.4, "n" => 6,
                              "shares" => Dict("0.5" => Dict("lo" => 0.25, "hi" => 0.45))),
                         Dict("fold_lo" => 0.5, "fold_hi" => 0.9, "n" => 6,
                              "shares" => Dict("0.5" => Dict("lo" => 0.55, "hi" => 0.75)))]
    end
    small = Dict("n" => 2,
                 "shares" => Dict("0.5" => Dict("lo" => 0.0, "hi" => 0.1)))
    t = Dict{String, Any}(
        "version" => "fixture", "spec" => "specs/155-fold-intervals",
        "catalog_build" => catalog_build, "solver_commit" => SERVED_SOLVER,
        "min_cell" => min_cell, "nominal_shares" => [0.5],
        "arms" => Dict("ax" => arm_b ? "b" : "a"),
        "table" => Dict("ax" => Dict(
            # the full-level cell is TOO SMALL, so a full-featured readout
            # must fall through to no_basis — the walk is what is tested
            "full" => Dict("down|False|direct|le3|True" => small),
            "no_basis" => Dict("down|False|le3|True" => cell),
            "no_converged" => Dict("down|False" => cell),
            "call_only" => Dict("down" => cell))))
    for k in drop
        delete!(t, k)
    end
    t
end

function write_table(t)
    path, io = mktemp()
    write(io, JSON3.write(t))
    close(io)
    path
end

load_fixture(t; served = SERVED_BUILD) =
    DeltaSignal.load_fold_ranges_table(write_table(t); served_build = served)

@testset "fold ranges (specs/155)" begin

    @testset "the loader's honest refusals" begin
        t, reason = DeltaSignal.load_fold_ranges_table(
            joinpath(mktempdir(), "nowhere.json"); served_build = SERVED_BUILD)
        @test t === nothing && reason == "missing"

        path, io = mktemp(); write(io, "{not json"); close(io)
        t, reason = DeltaSignal.load_fold_ranges_table(path; served_build = SERVED_BUILD)
        @test t === nothing && reason == "malformed"

        for key in ("table", "min_cell", "nominal_shares", "catalog_build", "arms")
            t, reason = load_fixture(fixture_table(drop = [key]))
            @test t === nothing && reason == "malformed"
        end

        t, reason = load_fixture(fixture_table(catalog_build = "someone-else"))
        @test t === nothing
        @test occursin("provenance mismatch", reason)
        @test occursin("someone-else", reason)

        for sentinel in (nothing, "unrecorded", "unreadable BUILD.json", "uninspectable")
            t, reason = load_fixture(fixture_table(); served = sentinel)
            @test t === nothing
            @test occursin("cannot be verified", reason)
        end

        t, reason = load_fixture(fixture_table())
        @test t !== nothing && reason == ""
    end

    @testset "the lookup walks the stratum hierarchy" begin
        t, _ = load_fixture(fixture_table())
        # fold 0.3: call down, |log2| ~1.74 -> le3; full cell too small, so
        # the no_basis cell serves, with ITS computed quantiles.
        r = DeltaSignal.readout_fold_ranges(t, "down", false, "direct", 0.3;
                                            converged = true)
        @test r !== nothing && length(r) == 1
        @test r[1]["axis"] == "ax" && r[1]["share"] == 0.5
        @test r[1]["lo"] == cell_lo && r[1]["hi"] == cell_hi
        @test r[1]["n"] == length(fixture_measured)

        # converged unknown: the converged-keyed levels are skipped and the
        # no_converged cell serves the same computed interval.
        r = DeltaSignal.readout_fold_ranges(t, "down", false, nothing, 0.3)
        @test r !== nothing && r[1]["lo"] == cell_lo

        # the endpoint's "normal" is the record's "not": no such cell here,
        # so nothing serves — never a borrowed cell.
        @test DeltaSignal.readout_fold_ranges(t, "normal", false, nothing, 0.98) === nothing

        # an untouched readout is served nothing (amendment-1 rule).
        @test DeltaSignal.readout_fold_ranges(t, "down", false, nothing, 1.0) === nothing

        # an up call has no cell at any level: honest absence.
        @test DeltaSignal.readout_fold_ranges(t, "up", false, nothing, 3.0) === nothing
    end

    @testset "arm (b): the band containing the fold refines; outside, the cell serves" begin
        t, _ = load_fixture(fixture_table(arm_b = true))
        inside = DeltaSignal.readout_fold_ranges(t, "down", false, nothing, 0.3;
                                                 converged = true)
        @test inside[1]["lo"] == 0.25 && inside[1]["hi"] == 0.45 && inside[1]["n"] == 6

        second = DeltaSignal.readout_fold_ranges(t, "down", false, nothing, 0.6;
                                                 converged = true)
        @test second[1]["lo"] == 0.55 && second[1]["hi"] == 0.75

        # 0.45 falls in the gap between bands: the cell's own interval, never
        # an extrapolated band.
        gap = DeltaSignal.readout_fold_ranges(t, "down", false, nothing, 0.45;
                                              converged = true)
        @test gap[1]["lo"] == cell_lo && gap[1]["hi"] == cell_hi
    end

    @testset "the shipped table is the record's, for the served catalog" begin
        shipped_path = joinpath(@__DIR__, "..", "src", "api", "fold_ranges.json")
        @test isfile(shipped_path)
        t, reason = DeltaSignal.load_fold_ranges_table(shipped_path;
                                                       served_build = SERVED_BUILD)
        @test t !== nothing && reason == ""
        @test string(t[:spec]) == "specs/155-fold-intervals"
        @test string(t[:solver_commit]) == SERVED_SOLVER
        @test Int(t[:min_cell]) == 30
        # the measured verdict: arm (a) on every axis — (b) narrowed nothing
        for (_, arm) in pairs(t[:arms])
            @test string(arm) == "a"
        end
        # every shipped cell honours the minimum and is a real interval
        for (_, levels) in pairs(t[:table]), (_, cells) in pairs(levels),
            (_, cell) in pairs(cells)
            @test Int(cell[:n]) >= 30
            for (_, iv) in pairs(cell[:shares])
                @test Float64(iv[:lo]) <= Float64(iv[:hi])
            end
        end
        # a different serving record refuses the shipped table too
        t2, reason2 = DeltaSignal.load_fold_ranges_table(shipped_path;
                                                         served_build = "other-build")
        @test t2 === nothing && occursin("provenance mismatch", reason2)
    end
end

println("test_fold_ranges.jl passed")
