#!/usr/bin/env julia
# Worked micro-example from the spec (Section 16) — math sanity check.
#
# Setup: one reaction with two activators A, B; one inhibitor I; one substrate S;
# product P (not used in SS). Sensitivity off (s=0). Hill h=2, K=0.05.
# Observed: A=0.5, B=0.3, I=0.2, S=0.2.
# Expected:
#   - Per-input ~x = x (s=0 ⇒ α=1)
#   - A_agg = geomean(0.5, 0.3) = √0.15 ≈ 0.387
#   - H = 1/(1 + 3·0.2²) = 1/1.12 ≈ 0.893
#   - L = S = 0.2 (single substrate; geomean of one input is itself)
#   - s_r = A·H·L ≈ 0.387·0.893·0.2 ≈ 0.069
#   - y = Hill(0.069, 2, 0.05) ≈ 0.655

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Test

push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
include("../src/core/sensitivity.jl")
include("../src/core/aggregators.jl")
include("../src/core/hill_functions.jl")

@testset "Spec Section 16 worked example" begin
    A_in = [0.5, 0.3]
    weights = [0.5, 0.5]
    A_agg = geometric_mean_aggregator(A_in, weights)
    @test isapprox(A_agg, sqrt(0.15); atol=1e-4)

    H = inhibition_aggregator([0.2], [3.0], [2.0])
    @test isapprox(H, 1.0 / (1.0 + 3.0 * 0.2^2); atol=1e-8)

    L = substrate_availability_aggregator([0.2], [1.0])
    @test isapprox(L, 0.2; atol=1e-4)

    s_r = A_agg * H * L
    y = hill_activation(s_r, 2.0, 0.05)

    println("worked-example pre-activation s_r = ", round(s_r; digits=4))
    println("worked-example output y           = ", round(y; digits=4))

    # Spec gives y ≈ 0.655. Allow some slack for the small ε log-safety term
    # in geomean which nudges A_agg slightly above pure sqrt(0.15).
    @test isapprox(s_r, 0.069; atol=2e-3)
    @test isapprox(y, 0.655; atol=5e-3)

    # Sensitivity transform with s=0 should be identity.
    for x in (0.1, 0.5, 0.9)
        @test apply_sensitivity_transform(x; s=0.0) == x
    end

    # Hill self-consistency at baseline x₀=0.01 with default h=2, K=0.1:
    # Hill(0.01, 2, 0.1) should be ≈ 0.01 so the baseline is a fixed point.
    @test isapprox(hill_activation(0.01, 2.0, 0.1), 0.01; atol=2e-3)
end
