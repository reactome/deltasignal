#!/usr/bin/env julia
#
# The API's error sanitizer, and the maintenance hazard inside it.
#
# `user_facing_error` deliberately does NOT echo arbitrary client input back
# into a response body, so it allow-lists the messages it is willing to
# forward and replaces everything else with a generic "Invalid input.".
#
# The hazard: adding a new `*_ERROR` constant and throwing it is not enough.
# If it is missing from the allow-list, the specific message is thrown,
# discarded, and silently downgraded. OBS_UNKNOWN_ERROR was stranded that way
# from the moment it was written. It names the catalog-mismatch failure mode
# -- observations whose uuids come from a different catalog build than the
# network -- which is the hardest failure in this project to diagnose from
# outside, and the client only ever saw "Invalid input.".
#
# The first test below is the guard: it enumerates the constants rather than
# listing them by hand, so a newly added one fails here until it is forwarded.

using Test

include(joinpath(@__DIR__, "..", "src", "DeltaSignal.jl"))
using .DeltaSignal

const SRV = DeltaSignal

@testset "API error sanitization" begin

    @testset "every error constant reaches the client" begin
        # Enumerated from the module, NOT hand-listed: a constant added later
        # is covered automatically.
        err_names = [n for n in Base.names(SRV; all=true)
                 if endswith(String(n), "_ERROR") && isdefined(SRV, n)]
        @test !isempty(err_names)

        for n in err_names
            msg = getfield(SRV, n)
            msg isa String || continue
            status, body = SRV.user_facing_error(ArgumentError(msg))
            @test status == 400
            @test body == msg
        end
    end

    @testset "the stranded one specifically" begin
        # Regression: this returned "Invalid input." and the catalog-mismatch
        # hint never reached anyone.
        _, body = SRV.user_facing_error(ArgumentError(SRV.OBS_UNKNOWN_ERROR))
        @test body == SRV.OBS_UNKNOWN_ERROR
        @test occursin("same catalog build", body)
    end

    @testset "arbitrary messages are still not echoed" begin
        # The reason the allow-list exists. Client-controlled text must not be
        # reflected into the response body.
        status, body = SRV.user_facing_error(ArgumentError("<script>alert(1)</script>"))
        @test status == 400
        @test body == "Invalid input."
        @test !occursin("script", body)

        _, body2 = SRV.user_facing_error(ArgumentError("pathway R-HSA-999999 is bad"))
        @test !occursin("R-HSA-999999", body2)
    end

    @testset "non-ArgumentError exceptions map to their own messages" begin
        @test SRV.user_facing_error(KeyError(:missing))[2] == "Required field missing from input."
        @test SRV.user_facing_error(SystemError("open"))[2] == "Could not read input file."
        for e in (KeyError(:x), SystemError("open"), ArgumentError("x"))
            @test SRV.user_facing_error(e)[1] == 400
        end
    end

    # The dev API served a two-month-old catalog while every benchmark ran
    # against newer ones, and nothing said so. /api/health now reports which
    # catalog is being served, and must not look healthy while serving nothing.
    @testset "health reports the catalog it is serving" begin
        mktempdir() do root
            # a recorded, complete build
            built = joinpath(root, "built")
            for p in ("R-HSA-1", "R-HSA-2")
                mkpath(joinpath(built, p))
                write(joinpath(built, p, "logic_network.csv"), "source_id,target_id\n")
            end
            write(joinpath(built, "BUILD.json"),
                  """{"build_id":"20260925-1200_abc1234","status":"complete",
                      "generator_commit":"abc1234","pathways_built":2}""")
            info = SRV.catalog_provenance(built)
            @test info["build_id"] == "20260925-1200_abc1234"
            @test info["generator_commit"] == "abc1234"
            @test info["status"] == "complete"
            @test info["pathways"] == 2
            @test !haskey(info, "warning")

            # a catalog with no manifest is reported as such, never guessed at
            legacy = joinpath(root, "legacy")
            mkpath(joinpath(legacy, "R-HSA-9"))
            write(joinpath(legacy, "R-HSA-9", "logic_network.csv"), "source_id,target_id\n")
            linfo = SRV.catalog_provenance(legacy)
            @test linfo["build_id"] == "unrecorded"
            @test linfo["pathways"] == 1

            # docker creates an EMPTY directory when a mount path is missing:
            # that must carry a warning, not read as a healthy catalog
            empty = joinpath(root, "empty"); mkpath(empty)
            einfo = SRV.catalog_provenance(empty)
            @test einfo["pathways"] == 0
            @test haskey(einfo, "warning")

            # a pathway directory without a network does not count as served
            partial = joinpath(root, "partial"); mkpath(joinpath(partial, "R-HSA-5"))
            @test SRV.catalog_provenance(partial)["pathways"] == 0

            # a path that does not exist at all
            @test SRV.catalog_provenance(joinpath(root, "nope"))["pathways"] == 0
        end
    end
end
