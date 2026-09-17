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
end
