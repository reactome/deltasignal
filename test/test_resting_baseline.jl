#!/usr/bin/env julia
#
# specs/124: DS_BASELINE, the internal level every parsed node rests at. The internal scale tops out at 1, so the
# ceiling on any node is 1 / DS_BASELINE fold. The default (0.01, a 100x ceiling) must leave every solve unchanged.

using Test

include(joinpath(@__DIR__, "..", "src", "DeltaSignal.jl"))
using .DeltaSignal

const EX = joinpath(@__DIR__, "..", "examples")

function with_baseline(f, value)
    prev = get(ENV, "DS_BASELINE", nothing)
    value === nothing ? delete!(ENV, "DS_BASELINE") : (ENV["DS_BASELINE"] = value)
    try
        return f()
    finally
        prev === nothing ? delete!(ENV, "DS_BASELINE") : (ENV["DS_BASELINE"] = prev)
    end
end

"""Two parents ANDed into T, built with the parser's resting level. Each parent is pinned at `fold` over its rest;
returns T's fold over rest and the solve's reported resting level."""
function and_two(fold::Float64)
    bl = DeltaSignal.resting_baseline()
    nodes = Dict(u => DeltaSignal.NetworkNode(u, u, "protein", nothing, u, bl) for u in ("A", "B", "T"))
    network = DeltaSignal.ReactionNetwork(nodes,
        [DeltaSignal.LogicNetworkEdge(p, "T", true, true, 1.0, "input") for p in ("A", "B")],
        Dict{String, DeltaSignal.SetExpansionMapping}())
    obs = Dict(p => (min(fold * bl, 1.0) * 100.0, 1.0) for p in ("A", "B"))   # UI = internal x 100
    params = DeltaSignal.SteadyStateParams(1.0, 0.1, 500, 1e-6, "penalty")
    result = DeltaSignal.solve_steady_state(network, obs, params)
    return result.node_activities["T"] / bl, result.diagnostics["resting_baseline"]
end

function sample_solve()
    network = DeltaSignal.parse_complete_network(joinpath(EX, "sample_logic_network.tsv"),
                                                 joinpath(EX, "sample_uuid_mapping.tsv"),
                                                 joinpath(EX, "sample_set_mappings.tsv"))
    first_node = sort(collect(keys(network.nodes)))[1]
    params = DeltaSignal.SteadyStateParams(1.0, 0.1, 500, 1e-6, "penalty")
    result = DeltaSignal.solve_steady_state(network, Dict(first_node => (0.0, 1.0)), params)
    return network, result
end

@testset "specs/124 resting baseline" begin
    @testset "default and validation" begin
        with_baseline(nothing) do
            @test DeltaSignal.resting_baseline() == 0.01
        end
        with_baseline("0.1") do
            @test DeltaSignal.resting_baseline() == 0.1
        end
        for bad in ("0", "-0.01", "0.6", "1", "abc", "")
            with_baseline(bad) do
                @test_throws ArgumentError DeltaSignal.resting_baseline()
            end
        end
    end

    @testset "the parser gives every node the resting level" begin
        for (value, expected) in ((nothing, 0.01), ("0.001", 0.001), ("0.1", 0.1))
            with_baseline(value) do
                network, _ = sample_solve()
                @test all(n.baseline == expected for n in values(network.nodes))
            end
        end
    end

    @testset "the default is today's solve, exactly" begin
        _, r0 = with_baseline(sample_solve, nothing)
        _, r1 = with_baseline(sample_solve, "0.01")
        @test r0.node_activities == r1.node_activities
        @test r0.diagnostics["resting_baseline"] == 0.01
    end

    @testset "the ceiling is 1 / baseline fold" begin
        # 5x AND 5x is 25x below every ceiling tested except 10x, where it stops at the top.
        with_baseline(nothing) do
            t, rb = and_two(5.0)
            @test isapprox(t, 25.0; rtol = 1e-6)
            @test rb == 0.01
        end
        with_baseline("0.001") do
            t, rb = and_two(5.0)
            @test isapprox(t, 25.0; rtol = 1e-6)
            @test rb == 0.001
        end
        with_baseline("0.1") do
            t, rb = and_two(5.0)
            @test t <= 10.0 + 1e-9
            @test isapprox(t, 10.0; rtol = 1e-3)
            @test rb == 0.1
        end
        # Near rest the multiplication does not depend on the ceiling.
        for value in (nothing, "0.1", "0.001")
            with_baseline(value) do
                t, _ = and_two(0.5)
                @test isapprox(t, 0.25; rtol = 1e-6)
            end
        end
    end
end
