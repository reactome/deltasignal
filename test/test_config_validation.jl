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
                 "DS_ASSEMBLY_LIMITING", "DS_OR_COMBINE", "DS_INHIBITOR_OR")
        haskey(ENV, name) && delete!(ENV, name)
    end
    config = resolve_reaction_eval_config()
    @test config.inhibition_mode == "divide"
    @test config.and_mode == "hill_log"
    @test config.or_mode == "mean"
    @test config.assembly_limiting == true
    @test config.or_combine == "max"
    @test config.inhibitor_or == false
    @test config.or_redundancy == 1.0
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
