#!/usr/bin/env julia
#
# The AND design intent, pinned at the magnitudes signals actually reach.
#
# Adam's specification: AND multiplies fold-changes with a Hill-like curve,
# constrained to 0-100 so that 100 x 100 gives 100 rather than 10,000, and
# "if it is near 1 I want it to be extremely close to pure multiplication" —
# 0.5*0.5 = 0.25, 1*1 = 1, 2*0.5 = 1.
#
# Every defect this feature fixes was EXACTLY IDENTITY at the 0.85/1.15
# classification cutoffs and wrong somewhere else in the range:
#
#   - `min` assembly clamp: identity for knockouts, discards every increase.
#   - `hill_log` at z_max=10: identity near baseline, reads a lone UI 50 as
#     41.4 and a lone UI 100 as 74.1, because e^10 lies far outside the UI
#     range so the tanh squashing is active throughout it.
#   - `hill_sat` at eps=1e-3: identity mid-range, but the epsilon exceeds the
#     internal values where knockouts live (~0.0006) and acts as a floor,
#     lifting 0.25*0.25 to 0.09 and biasing the whole model upward.
#
# Every check anyone ran was near the cutoffs, so all three survived. These
# assertions test where the failures live.

using Test

include(joinpath(@__DIR__, "..", "src", "DeltaSignal.jl"))
using .DeltaSignal
using .DeltaSignal: resolve_reaction_eval_config

"""Evaluate the configured AND mode on a set of fold-changes, in fold units."""
function and_fold(folds::Vector{Float64}; mode::String, eps::String="1e-5")
    prev_mode = get(ENV, "DS_AND_MODE", nothing)
    prev_eps = get(ENV, "DS_HILL_SAT_EPS", nothing)
    ENV["DS_AND_MODE"] = mode
    ENV["DS_HILL_SAT_EPS"] = eps
    try
        bl = 0.01
        network = DeltaSignal.ReactionNetwork(
            Dict(
                "T" => DeltaSignal.NetworkNode("T", "T", "protein", nothing, "T", bl),
                (("P$i" => DeltaSignal.NetworkNode("P$i", "P$i", "protein", nothing, "P$i", bl))
                 for i in eachindex(folds))...,
            ),
            [DeltaSignal.LogicNetworkEdge("P$i", "T", true, true, 1.0, "input")
             for i in eachindex(folds)],
            Dict{String, DeltaSignal.SetExpansionMapping}(),
        )
        # Pin each parent at its fold (UI = fold, since baseline UI is 1).
        obs = Dict("P$i" => (folds[i], 1.0) for i in eachindex(folds))
        params = DeltaSignal.SteadyStateParams(1.0, 0.1, 500, 1e-6, "penalty")
        result = DeltaSignal.solve_steady_state(network, obs, params)
        return result.node_activities["T"] / bl
    finally
        prev_mode === nothing ? delete!(ENV, "DS_AND_MODE") : (ENV["DS_AND_MODE"] = prev_mode)
        prev_eps === nothing ? delete!(ENV, "DS_HILL_SAT_EPS") : (ENV["DS_HILL_SAT_EPS"] = prev_eps)
    end
end

# A top-level `@testset` THROWS when it finishes with any failure, which
# aborts the file and silently skips every testset below it. That is how a
# config mismatch in the dev container hid 61 assertions here -- the cofactor
# tests, the silo-bridge tests and the observation-membership test all looked
# green because they never ran at all.
#
# Nesting them inside one outer testset makes failures accumulate: every
# testset executes, and the outer one reports the full tally at the end. A
# failing run still exits non-zero, it just stops lying about coverage.
@testset "AND curve behaviour" begin

@testset "AND is multiplication near baseline" begin
    # "extremely close to pure multiplication" where signals usually sit.
    @test and_fold([0.5, 0.5]; mode="hill_sat") ≈ 0.25 atol=0.01
    @test and_fold([1.0, 1.0]; mode="hill_sat") ≈ 1.00 atol=0.01
    @test and_fold([2.0, 0.5]; mode="hill_sat") ≈ 1.00 atol=0.01
    @test and_fold([2.0, 2.0]; mode="hill_sat") ≈ 4.00 atol=0.01
end

@testset "AND is constrained to the 0-100 range, not clipped early" begin
    # A value already INSIDE the range must pass through unchanged. This is
    # the requirement `hill_log` fails.
    @test and_fold([50.0]; mode="hill_sat") ≈ 50.0 atol=0.1
    @test and_fold([100.0]; mode="hill_sat") ≈ 100.0 atol=0.1
    # 10 x 10 = 100, not 10,000: saturation at the ceiling, not before it.
    @test and_fold([10.0, 10.0]; mode="hill_sat") ≈ 100.0 atol=0.1
    @test and_fold([100.0, 100.0]; mode="hill_sat") ≈ 100.0 atol=0.1
end

@testset "hill_log compresses inside the range (documents why the default changes)" begin
    # Not a defence of these values — a record of the behaviour being replaced,
    # so a future reader can see the difference rather than rediscover it.
    @test and_fold([50.0]; mode="hill_log") ≈ 41.4 atol=1.0
    @test and_fold([100.0]; mode="hill_log") ≈ 74.1 atol=1.0
    # It is identity near baseline, which is why every previous check passed.
    @test and_fold([0.5, 0.5]; mode="hill_log") ≈ 0.25 atol=0.01
end

@testset "hill_sat epsilon must not floor the low end" begin
    # eps=1e-3 exceeds the internal values where knockouts live (~0.0006), so
    # it lifts them and biases the model upward: measured, it predicted UP on
    # 317 of 564 cases against 246 actually UP.
    @test and_fold([0.25, 0.25]; mode="hill_sat", eps="1e-3") > 0.08
    @test and_fold([0.25, 0.25]; mode="hill_sat", eps="1e-5") ≈ 0.0625 atol=0.001
    # A knockout must still collapse the target.
    @test and_fold([0.0, 1.0]; mode="hill_sat", eps="1e-5") < 0.01
end

@testset "defaults are what measured best, which is NOT the design intent" begin
    # An honest name, because these two disagree and the tension is real.
    #
    # specs/002 R5 argued `hill_sat` implements Adam's stated AND intent --
    # "if it is near 1 I want it to be extremely close to pure multiplication"
    # -- and it does: 10x10 reads 99.98 under hill_sat against 74.06 under
    # hill_log, and a lone node at UI 50 reads 50.00 against 41.43. By that
    # argument hill_log is wrong, and the testsets above still pin those
    # curve shapes because the argument has not been withdrawn.
    #
    # But specs/002 chose on +31 over 564 experimental cases. Re-measured
    # 2026-09-16 on 23,022 wide-curator cases, one catalog build, one variable
    # at a time: `hill_sat` costs -90 and `assembly_limiting=false` costs -196.
    # FR-008 makes the wide set the decision basis, so the defaults follow the
    # measurement and this test records that they are not the design intent.
    #
    # What that means: the compression hill_log applies inside the operating
    # range is apparently doing useful work that "correct" multiplication does
    # not. Nobody has explained why. Until someone does, this is an empirical
    # default, not a principled one -- see specs/009-solver-defaults.
    for name in ("DS_AND_MODE", "DS_HILL_SAT_EPS", "DS_ASSEMBLY_LIMITING")
        haskey(ENV, name) && delete!(ENV, name)
    end
    config = resolve_reaction_eval_config()
    @test config.and_mode == "hill_log"
    @test config.assembly_limiting == true
    # Inert while and_mode is hill_log; kept correctly sized so switching the
    # AND mode cannot silently restore a 10%-of-baseline epsilon.
    @test config.hill_sat_eps ≈ 1e-5
end


end  # outer testset
