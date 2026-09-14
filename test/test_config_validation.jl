#!/usr/bin/env julia
#
# Regression tests for configuration guard rails and the API's silent
# network-substitution fallbacks.
#
# Uses `Test` and real `@test` assertions, so a regression makes this file exit
# non-zero. Most other files in this directory wrap everything in try/catch,
# print a failure and return `false` to a caller that discards it, so they
# cannot fail — see test/README or the suite notes.

using Test

include(joinpath(@__DIR__, "..", "src", "DeltaSignal.jl"))
using .DeltaSignal
using .DeltaSignal: resolve_reaction_eval_config, DS_VALID_MODES,
                    _bool_env, _float_env, _mode_env

"""Run `f` with `ENV[name] = value`, restoring the previous state after."""
function with_env(f, name::String, value)
    had = haskey(ENV, name)
    old = had ? ENV[name] : ""
    value === nothing ? delete!(ENV, name) : (ENV[name] = value)
    try
        f()
    finally
        had ? (ENV[name] = old) : delete!(ENV, name)
    end
end

@testset "DS_* mode validation" begin
    # Every value the dispatch actually implements must still be accepted.
    # This is the half that matters: an over-tight allowlist breaks working
    # configurations, which is worse than the bug it replaces.
    @testset "all implemented modes accepted" begin
        for (name, values) in DS_VALID_MODES, value in values
            with_env(name, value) do
                @test resolve_reaction_eval_config() isa DeltaSignal.ReactionEvalConfig
            end
        end
    end

    # An unrecognised mode used to fall through to a *different model*. For
    # DS_INHIBITION_MODE the fallback was "spec" with default beta 0, i.e.
    # inhibition silently switched off entirely.
    @testset "typos rejected, not silently remapped" begin
        for bad in ("Divide", "divide ", "dividee", "", "DIVIDE")
            with_env("DS_INHIBITION_MODE", bad) do
                @test_throws ArgumentError resolve_reaction_eval_config()
            end
        end
        for bad in ("Hill_Log", "hill-log", "bogus")
            with_env("DS_AND_MODE", bad) do
                @test_throws ArgumentError resolve_reaction_eval_config()
            end
        end
        with_env("DS_OR_COMBINE", "Gate") do
            @test_throws ArgumentError resolve_reaction_eval_config()
        end
    end

    @testset "error names the variable and lists valid values" begin
        with_env("DS_INHIBITION_MODE", "nope") do
            err = try resolve_reaction_eval_config() catch e; e end
            @test err isa ArgumentError
            @test occursin("DS_INHIBITION_MODE", err.msg)
            @test occursin("divide", err.msg)
        end
    end
end

@testset "DS_* boolean parsing" begin
    # Flags were tested two incompatible ways: some `== "1"`, others `!= "0"`.
    # DS_ASSEMBLY_LIMITING="true" therefore turned the feature OFF while
    # DS_SCC_SOLVE="false" left it ON. Both spellings must now mean the same.
    @testset "true spellings" begin
        for value in ("1", "true", "TRUE", "yes", "on", " true ")
            @test _bool_env("DS_TEST_FLAG", false) == false  # unset -> default
            with_env("DS_TEST_FLAG", value) do
                @test _bool_env("DS_TEST_FLAG", false) == true
            end
        end
    end

    @testset "false spellings" begin
        for value in ("0", "false", "FALSE", "no", "off")
            with_env("DS_TEST_FLAG", value) do
                @test _bool_env("DS_TEST_FLAG", true) == false
            end
        end
    end

    @testset "non-boolean rejected" begin
        for value in ("bogus", "2", "")
            with_env("DS_TEST_FLAG", value) do
                @test_throws ArgumentError _bool_env("DS_TEST_FLAG", false)
            end
        end
    end

    @testset "assembly limiting honours word spellings" begin
        with_env("DS_ASSEMBLY_LIMITING", "true") do
            @test resolve_reaction_eval_config().assembly_limiting == true
        end
        with_env("DS_ASSEMBLY_LIMITING", "false") do
            @test resolve_reaction_eval_config().assembly_limiting == false
        end
    end
end

@testset "DS_* numeric parsing" begin
    @test _float_env("DS_TEST_NUM", 2.5) == 2.5           # unset -> default
    with_env("DS_TEST_NUM", "3.5") do
        @test _float_env("DS_TEST_NUM", 2.5) == 3.5
    end
    for bad in ("abc", "", "NaN", "Inf")
        with_env("DS_TEST_NUM", bad) do
            @test_throws ArgumentError _float_env("DS_TEST_NUM", 2.5)
        end
    end
    with_env("DS_HILL_LOG_ZMAX", "not-a-number") do
        err = try resolve_reaction_eval_config() catch e; e end
        @test err isa ArgumentError
        @test occursin("DS_HILL_LOG_ZMAX", err.msg)
    end
end

@testset "defaults are the validated winning config" begin
    for name in ("DS_INHIBITION_MODE", "DS_AND_MODE", "DS_OR_MODE",
                 "DS_ASSEMBLY_LIMITING", "DS_OR_COMBINE", "DS_INHIBITOR_OR",
                 "DS_HILL_SAT_EPS")
        haskey(ENV, name) && delete!(ENV, name)
    end
    config = resolve_reaction_eval_config()
    @test config.inhibition_mode == "divide"
    # hill_sat / clamp-off / eps 1e-5 as of specs/002-upregulation-propagation.
    # test/test_and_curves.jl covers why; this only pins that the defaults are
    # what that feature measured.
    @test config.and_mode == "hill_sat"
    @test config.or_mode == "mean"
    @test config.assembly_limiting == false
    @test config.hill_sat_eps == 1e-5
    # A divide-by-zero guard, and nothing else. Baseline is 0.01, so this
    # must stay orders of magnitude below it — at the old 1e-3 it was 10% of
    # baseline and silently set the de-repression ceiling, compressed the
    # response curve and shifted maximum suppression.
    @test config.inhibitor_eps == 1e-12
    @test config.inhibitor_eps < 0.01 / 1000
    @test config.or_combine == "max"
    @test config.inhibitor_or == false
    @test config.or_redundancy == 1.0
end

@testset "the inhibition epsilon is a guard, not a model parameter" begin
    # The defect this pins: at 1e-3 the epsilon was 10% of baseline and did
    # three things nobody chose. At a true guard size it does one.
    bl = 0.01
    network = DeltaSignal.ReactionNetwork(
        Dict(
            "T" => DeltaSignal.NetworkNode("T", "T", "protein", nothing, "T", bl),
            "I" => DeltaSignal.NetworkNode("I", "I", "protein", nothing, "I", bl),
        ),
        [DeltaSignal.LogicNetworkEdge("I", "T", true, false, 1.0, "regulator")],
        Dict{String, DeltaSignal.SetExpansionMapping}(),
    )
    function target_fold(eps, inhibitor_ui)
        prev = get(ENV, "DS_INHIBITOR_EPS", nothing)
        ENV["DS_INHIBITOR_EPS"] = string(eps)
        try
            params = DeltaSignal.SteadyStateParams(1.0, 0.1, 500, 1e-6, "penalty")
            result = DeltaSignal.solve_steady_state(network, Dict("I" => (inhibitor_ui, 1.0)), params)
            return result.node_activities["T"] / bl
        finally
            prev === nothing ? delete!(ENV, "DS_INHIBITOR_EPS") : (ENV["DS_INHIBITOR_EPS"] = prev)
        end
    end

    # Halving the inhibitor should double its target: h = bl/x exactly.
    # The old default returned 1.833 here — an 8% error from a "guard".
    @test target_fold(1e-12, 0.5) ≈ 2.0 atol=0.01
    @test target_fold(1e-3, 0.5) < 1.9        # the defect, pinned

    # A full knockout is bounded by the per-reaction clamp, not by epsilon,
    # so shrinking epsilon by nine orders of magnitude must not change it.
    @test target_fold(1e-12, 0.0) ≈ target_fold(1e-9, 0.0) atol=1e-6
    @test target_fold(1e-12, 0.0) <= 10.0 + 1e-6
end

@testset "damping is range-checked" begin
    # nv = (1-λ)x + λF(x) is written back unclamped, so λ outside [0,1] left the
    # model's [0,1] domain and threw an uncaught DomainError from inside
    # hill_log; λ = 0 silently never updated the component.
    network = DeltaSignal.ReactionNetwork(
        Dict(
            "A" => DeltaSignal.NetworkNode("A", "A", "protein", nothing, "A", 0.01),
            "B" => DeltaSignal.NetworkNode("B", "B", "protein", nothing, "B", 0.01),
        ),
        [DeltaSignal.LogicNetworkEdge("A", "B", true, true, 1.0, "input")],
        Dict{String, DeltaSignal.SetExpansionMapping}(),
    )
    observations = Dict("A" => (50.0, 1.0))
    params = DeltaSignal.SteadyStateParams(1.0, 0.1, 100, 1e-6, "penalty")

    with_env("DS_SCC_DAMPING", "0.5") do
        @test DeltaSignal.solve_steady_state(network, observations, params) !== nothing
    end
    for bad in ("1.5", "-0.5", "0.0", "3.0")
        with_env("DS_SCC_DAMPING", bad) do
            @test_throws ArgumentError DeltaSignal.solve_steady_state(
                network, observations, params)
        end
    end
    for bad in ("1.0", "-0.1")
        with_env("DS_DAMPING", bad) do
            @test_throws ArgumentError DeltaSignal.solve_steady_state(
                network, observations, params)
        end
    end
end
