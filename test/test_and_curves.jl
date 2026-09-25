#!/usr/bin/env julia
#
# The AND design intent, pinned at the magnitudes signals actually reach.
#
# Specification: AND multiplies fold-changes with a Hill-like curve,
# constrained to 0-100 so that 100 x 100 gives 100 rather than 10,000, and
# near baseline it must behave as close to pure multiplication as possible —
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

"""Evaluate an AND mode on fold-changes. `eps` defaults to the SHIPPED
DS_HILL_SAT_EPS so callers that omit it exercise real behaviour; it was
pinned at "1e-5" after the code default moved, which made every test that
omitted it measure a configuration the code no longer ships."""
function and_fold(folds::Vector{Float64}; mode::String, eps::String="1e-9")
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
    # As close to pure multiplication as possible where signals usually sit.
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

@testset "defaults ARE the design intent (specs/010)" begin
    # An honest name, because these two disagree and the tension is real.
    #
    # specs/002 R5 argued `hill_sat` implements the stated AND intent -- near
    # baseline it should behave as close to pure multiplication as possible
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
    @test config.and_mode == "hill_sat"
    # Off since specs/023: under root pinning min() blocked every single-subunit
    # overexpression (held-out +401 with it off).
    @test config.assembly_limiting == false
    # Inert while and_mode is hill_log; kept correctly sized so switching the
    # AND mode cannot silently restore a 10%-of-baseline epsilon.
    @test config.hill_sat_eps ≈ 1e-9
end

@testset "hill_log_asym multiplies faithfully below baseline" begin
    # AND is specified as multiplication of fold-changes. hill_log applies
    # z_max*tanh(log_fold/z_max) SYMMETRICALLY, so it compresses downward as
    # hard as upward -- but there is nothing to saturate against downward.
    # The internal domain is [0,1] with baseline 0.01, so an upward fold is
    # genuinely capped at 100x (a real ceiling, which tanh models), while a
    # downward fold of 0.01 is perfectly representable. Compressing it is not
    # modelling a floor, it is error.
    #
    # Measured error of hill_log against the product it should compute:
    #   0.5 x 2.0 -> 0.0%   |  0.25 x 0.25 ->  +7.2%
    #   2.0 x 2.0 -> -0.9%  |  0.1  x 0.1  -> +35.2%
    #   0.5 x 0.5 -> +0.9%  |  0.001 x 1.0 -> +167.6%
    # Accurate above baseline, unbounded and systematically UPWARD below it,
    # which biases every down-regulated value toward baseline.

    asym(folds) = and_fold(folds; mode="hill_log_asym")

    @testset "sub-baseline products are exact" begin
        @test asym([0.5, 0.5]) ≈ 0.25 rtol=1e-6
        @test asym([0.25, 0.25]) ≈ 0.0625 rtol=1e-6
        @test asym([0.1, 0.1]) ≈ 0.01 rtol=1e-6
        @test asym([0.5, 0.5, 0.5]) ≈ 0.125 rtol=1e-6
        @test asym([0.01, 1.0]) ≈ 0.01 rtol=1e-6
        # Two orders below baseline, where hill_log reads +167.6% high.
        @test asym([0.001, 1.0]) ≈ 0.001 rtol=1e-6
    end

    @testset "a zero input gives exactly zero" begin
        # hill_log's hardcoded eps=1e-6 against baseline 0.01 made a knockout
        # contribute log(1e-6/0.010001) = -9.21 instead of -Inf, so 0 x
        # anything came out at 0.0007 and rose with the co-input. Same epsilon
        # sizing bug already fixed in DS_HILL_SAT_EPS and DS_INHIBITOR_EPS.
        @test asym([0.0, 1.0]) == 0.0
        @test asym([0.0, 100.0]) == 0.0
        @test asym([0.0, 0.0]) == 0.0
        # A required input at zero cannot be rescued by ANY co-input.
        for co in (1.0, 10.0, 100.0)
            @test asym([0.0, co]) == 0.0
        end
    end

    @testset "crossing baseline is exact" begin
        @test asym([0.5, 2.0]) ≈ 1.0 rtol=1e-6
        @test asym([0.1, 10.0]) ≈ 1.0 rtol=1e-6
        @test asym([1.0, 1.0]) ≈ 1.0 rtol=1e-9
    end

    @testset "upward compression is retained" begin
        # The point of the asymmetry: the 100x ceiling is real, so the upward
        # branch keeps hill_log's behaviour rather than reverting to raw
        # multiplication. These are hill_log's numbers, not the products.
        @test asym([2.0, 2.0]) ≈ 3.9649 rtol=1e-3      # not 4.0
        @test asym([10.0, 10.0]) ≈ 74.07 rtol=1e-3     # not 100
        @test asym([100.0, 100.0]) ≈ 100.0 rtol=1e-6   # capped, not 10000
        @test asym([10.0, 10.0]) < 100.0
    end

    @testset "strictly better than hill_log below baseline" begin
        # Direct comparison on the same inputs: the asymmetric mode is closer
        # to the product everywhere below baseline.
        for folds in ([0.5, 0.5], [0.25, 0.25], [0.1, 0.1], [0.001, 1.0])
            want = prod(folds)
            @test abs(asym(folds) - want) <= abs(and_fold(folds; mode="hill_log") - want)
        end
    end
end

@testset "hill_sat implements the stated AND intent exactly" begin
    # The specification: AND multiplies fold-changes, capped at 100, and
    # `1/2 * 1/2` should be close to `1/4` -- below baseline matters as much as
    # above it. hill_sat is the mode that implements that shape; it had two
    # epsilons standing in the way, both now removed.
    sat(folds) = and_fold(folds; mode="hill_sat")

    @testset "products are exact below baseline" begin
        @test sat([0.5, 0.5]) ≈ 0.25 rtol=1e-4
        @test sat([0.25, 0.25]) ≈ 0.0625 rtol=1e-4
        @test sat([0.1, 0.1]) ≈ 0.01 rtol=1e-4
        @test sat([0.5, 0.5, 0.5]) ≈ 0.125 rtol=1e-4
        # Where the old floors bit hardest: eps=1e-5 read these 20.7% and
        # 452% high respectively.
        @test sat([0.001, 1.0]) ≈ 0.001 rtol=1e-3
        @test sat([0.0001, 1.0]) ≈ 0.0001 rtol=1e-3
    end

    @testset "a zero input gives exactly zero" begin
        # The ratio carried a hardcoded eps of 1e-6 against baseline 0.01, so
        # a knockout contributed a fold of 1e-4 rather than 0 and `0 x 100`
        # came out at 0.0100 -- rising with the co-input, so an abundant
        # partner could "rescue" a knockout.
        @test sat([0.0, 1.0]) == 0.0
        @test sat([0.0, 100.0]) == 0.0
        @test sat([0.0, 0.0]) == 0.0
    end

    @testset "the 100x ceiling is real and still applies" begin
        @test sat([10.0, 10.0]) ≈ 100.0 rtol=1e-3     # product is exactly 100
        @test sat([100.0, 100.0]) ≈ 100.0 rtol=1e-6   # 10000 capped to 100
        @test sat([2.0, 2.0]) ≈ 4.0 rtol=1e-3
        @test sat([0.5, 2.0]) ≈ 1.0 rtol=1e-6
    end

    @testset "closer to the product than hill_log everywhere" begin
        for folds in ([0.5,0.5], [0.1,0.1], [0.001,1.0], [10.0,10.0], [2.0,2.0])
            want = min(prod(folds), 100.0)
            @test abs(sat(folds) - want) <= abs(and_fold(folds; mode="hill_log") - want)
        end
    end
end

@testset "wide AND reactions saturate instead of collapsing" begin
    # A numerical-stability bug, not a modelling one, and severe: the smooth-min
    # capping hill_sat at 100x was written as
    #     (raw + max - sqrt(d^2 + eps^2)) / 2,  d = raw - max
    # which CATASTROPHICALLY CANCELS once raw is large. With 100 AND inputs at
    # fold 2, raw is 1.3e28; sqrt(d^2) equals d to machine precision, so
    # (raw + 1 - (raw - 1))/2 evaluates to 0 rather than 1, and the smooth-max
    # below then returned eps/2. A strongly ELEVATED wide reaction therefore
    # read as ~0 -- the saturation inverted.
    #
    # Class_I_MHC has reactions carrying hundreds of nodes and did not solve at
    # all under hill_sat because of this. No assertion covered wide reactions,
    # so it survived; these are the ones that would have caught it.
    wide(n, f) = and_fold(fill(f, n); mode="hill_sat")

    @testset "elevated wide reactions cap at 100x" begin
        for n in (10, 50, 100, 200, 400)
            @test wide(n, 2.0) ≈ 100.0 rtol=1e-6
        end
        # The failure signature was collapse to the epsilon floor, orders
        # BELOW baseline, for something that should read maximally elevated.
        @test wide(400, 2.0) > 1.0
        @test wide(400, 10.0) ≈ 100.0 rtol=1e-6
    end

    @testset "wide reactions at baseline stay at baseline" begin
        # A product of many 1.0s is 1.0; drift here would move every node in a
        # hub reaction.
        for n in (10, 100, 400)
            @test wide(n, 1.0) ≈ 1.0 rtol=1e-9
        end
    end

    @testset "suppressed wide reactions stay far below baseline" begin
        # These floor at the smooth-max epsilon rather than reaching the true
        # product (0.5^50 = 8.9e-16). That floor is ~5e-8 in fold terms, seven
        # orders below the 0.85 DOWN cutoff, so it cannot change a call.
        for n in (10, 50, 200)
            @test wide(n, 0.5) < 0.01
        end
    end
end

end  # outer testset
