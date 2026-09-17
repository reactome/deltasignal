#!/usr/bin/env julia
#
# Observation pinning semantics, and the confidence gate in particular.
#
# `solve_steady_state` takes observations as (activity, confidence) and the
# pinning step gates on confidence: anything at or below OBS_CONFIDENCE_TOL is
# meant to be "no observation at all". That gate was inoperative for ROOT
# nodes, and roots are the case that matters — a perturbed gene usually has no
# incoming reaction.
#
# The mechanism: the initial guess seeded x0 from the observation without
# applying the gate. A root has no incoming reaction, so nothing ever
# overwrote that seed and a confidence-0 observation was enforced exactly as
# hard as a confidence-1 one. Non-root nodes hid it completely, because the
# forward model overwrites them on the first sweep — which is why every check
# anyone ran looked right.
#
# The benchmark always sends confidence 1.0, so no measured number moved; this
# is a correctness defect on the documented input contract, which ships
# confidences of 0.9/0.8/0.7 in examples/sample_observations.csv.

using Test

include(joinpath(@__DIR__, "..", "src", "DeltaSignal.jl"))
using .DeltaSignal
using .DeltaSignal: OBS_CONFIDENCE_TOL

const BL = 0.01
const P = DeltaSignal.SteadyStateParams(1.0, 0.1, 500, 1e-6, "penalty")

mknet(nodes, edges) = DeltaSignal.ReactionNetwork(
    Dict(n => DeltaSignal.NetworkNode(n, n, "protein", nothing, n, BL) for n in nodes),
    edges, Dict{String, DeltaSignal.SetExpansionMapping}())

edge(s, t) = DeltaSignal.LogicNetworkEdge(s, t, true, true, 1.0, "input")

"""Solve and return the node's value as a fold of baseline."""
function fold(net, obs, node)
    r = DeltaSignal.solve_steady_state(net, obs, P)
    return r.node_activities[node] / BL
end

# One outer testset so a failure in an early group does not abort the file and
# silently skip every group below it.
@testset "observation pinning" begin

    @testset "confidence gate on a ROOT node" begin
        # S has no incoming reaction — the case the gate failed on.
        net = mknet(["S", "T"], [edge("S", "T")])

        # Confident observations pin, and the pinned value is exact.
        @test fold(net, Dict("S" => (100.0, 1.0)), "S") ≈ 100.0
        @test fold(net, Dict("S" => (100.0, 0.5)), "S") ≈ 100.0
        @test fold(net, Dict("S" => (100.0, 1.0)), "T") > 1.0

        # At or below the tolerance the observation must not be applied at all.
        for conf in (0.0, 1e-9, OBS_CONFIDENCE_TOL)
            @test fold(net, Dict("S" => (100.0, conf)), "S") ≈ 1.0
            @test fold(net, Dict("S" => (100.0, conf)), "T") ≈ 1.0
        end

        # A knockout must be dropped by the gate too, not just an increase —
        # otherwise the bug survives in the direction the benchmark uses most.
        @test fold(net, Dict("S" => (0.0, 0.0)), "S") ≈ 1.0
        @test fold(net, Dict("S" => (0.0, 1.0)), "S") ≈ 0.0
    end

    @testset "confidence gate on a NON-ROOT node" begin
        # M is downstream of U, so the forward model overwrites it. This path
        # always worked; it is here so a future change cannot fix roots by
        # breaking this.
        net = mknet(["U", "M", "T"], [edge("U", "M"), edge("M", "T")])
        @test fold(net, Dict("M" => (100.0, 1.0)), "M") ≈ 100.0
        @test fold(net, Dict("M" => (100.0, 0.0)), "M") ≈ 1.0
    end

    @testset "the gate is one threshold, not two" begin
        # The initial guess and the pinning step each had their own literal and
        # disagreed. Both now read OBS_CONFIDENCE_TOL; a value just above it
        # must pin, and the boundary itself must not.
        net = mknet(["S", "T"], [edge("S", "T")])
        @test fold(net, Dict("S" => (100.0, OBS_CONFIDENCE_TOL)), "S") ≈ 1.0
        @test fold(net, Dict("S" => (100.0, nextfloat(OBS_CONFIDENCE_TOL))), "S") ≈ 100.0
    end

    @testset "unobserved nodes stay exactly at baseline" begin
        # A node off the perturbed path must not drift. Drift here would read
        # as a false change, which is the single largest error category.
        net = mknet(["S", "X", "U", "V"], [edge("S", "X"), edge("U", "V")])
        @test fold(net, Dict("S" => (100.0, 1.0)), "V") ≈ 1.0
        @test fold(net, Dict("S" => (100.0, 1.0)), "U") ≈ 1.0
    end

    @testset "an all-baseline reaction returns baseline at any width" begin
        # If a reaction whose inputs are all at baseline did not return
        # baseline, every unperturbed node downstream of a hub would drift.
        for n in (1, 2, 3, 5, 8, 20)
            parents = ["P$i" for i in 1:n]
            net = mknet(vcat(parents, ["T"]), [edge(p, "T") for p in parents])
            obs = Dict(p => (1.0, 1.0) for p in parents)
            @test fold(net, obs, "T") ≈ 1.0
        end
    end
end
