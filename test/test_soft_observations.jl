#!/usr/bin/env julia
#
# specs/156: soft observations — a node whose solve target becomes
# `w*y + (1-w)*F(x)`. The fixtures COMPUTE the blend's fixed points from
# networks whose unobserved answers the test also computes (never pinning the
# quantity the rule reasons about): weight 1 equals the hard pin, weight at
# the gate is discarded, in between the value is monotone in w and strictly
# between the two ends, downstream reads the blend, roots settle at theirs,
# precedence is reported, and a request without soft observations is
# BIT-IDENTICAL whatever DS_SOFT_OBS says.

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

solve(net, obs; soft = Dict{String, Tuple{Float64, Float64}}()) =
    DeltaSignal.solve_steady_state(net, obs, P; soft_observations = soft)

act(r, n) = r.node_activities[n]

@testset "soft observations (specs/156)" begin

    # A → B → C. A pinned up; B soft-observed BELOW what the network computes.
    net = mknet(["A", "B", "C"], [edge("A", "B"), edge("B", "C")])
    obs = Dict("A" => (50.0, 1.0))
    soft_val = 2.0   # UI scale: 2x baseline, far below what A=50x drives

    @testset "weight 1 equals the hard pin; the gate discards" begin
        free = solve(net, obs)
        hard = solve(net, merge(obs, Dict("B" => (soft_val, 1.0))))
        softened = solve(net, obs; soft = Dict("B" => (soft_val, 1.0)))
        @test act(softened, "B") ≈ act(hard, "B") atol = 1e-9
        @test act(softened, "C") ≈ act(hard, "C") rtol = 1e-6
        @test softened.converged

        for w in (0.0, 1e-9, OBS_CONFIDENCE_TOL)
            gated = solve(net, obs; soft = Dict("B" => (soft_val, w)))
            @test act(gated, "B") ≈ act(free, "B") rtol = 1e-9
            @test gated.diagnostics["soft_observed"] == 0
            @test gated.diagnostics["soft_obs_rule"] == "off"
        end
    end

    @testset "between the ends: strictly between, monotone in w, downstream reads it" begin
        free = solve(net, obs)
        hard_b = soft_val / 100.0
        free_b = act(free, "B")
        @test free_b > hard_b   # the fixture computes: the network disagrees upward

        last_b = free_b
        last_c = act(free, "C")
        for w in (0.25, 0.5, 0.75)
            r = solve(net, obs; soft = Dict("B" => (soft_val, w)))
            b = act(r, "B")
            # strictly between the unobserved answer and the observation
            @test hard_b < b < free_b
            # monotone: more weight pulls B further toward the observation
            @test b < last_b
            # downstream C reads the reconciled value, so it falls too
            @test act(r, "C") < last_c
            @test r.diagnostics["soft_observed"] == 1
            @test r.diagnostics["soft_obs_rule"] == "blend"
            # the blend's own fixed point: B == w*y + (1-w)*F(B's inputs),
            # where F here is what A alone drives — computed, not assumed:
            # with A pinned, F(B) is the free solve's B.
            @test b ≈ w * hard_b + (1 - w) * free_b rtol = 1e-6
            last_b, last_c = b, act(r, "C")
        end
    end

    @testset "a soft-observed ROOT settles at its blend with baseline" begin
        # S has no producing reaction; its F is its resting value (1x).
        net2 = mknet(["S", "T"], [edge("S", "T")])
        for w in (0.25, 0.75)
            r = solve(net2, Dict{String, Tuple{Float64, Float64}}();
                      soft = Dict("S" => (100.0, w)))
            @test act(r, "S") / BL ≈ w * 100.0 + (1 - w) * 1.0 rtol = 1e-6
            @test act(r, "T") / BL > 1.0   # downstream reads the lifted root
            # the seed IS the root's fixed point: the solve must say so (the
            # final check used to re-blend the root's own value and report a
            # phantom residual w*(1-w)*|y - x| — the first arm run's killer)
            @test r.converged
            @test r.diagnostics["soft_observed"] == 1
        end
    end

    @testset "precedence: a hard pin wins and is reported" begin
        r = solve(net, merge(obs, Dict("B" => (80.0, 1.0)));
                  soft = Dict("B" => (soft_val, 0.5)))
        @test act(r, "B") ≈ 0.80 atol = 1e-9
        @test r.diagnostics["soft_observed"] == 0
        over = r.diagnostics["soft_overridden"]
        @test length(over) == 1
        @test over[1]["uuid"] == "B" && over[1]["reason"] == "pin"
    end

    @testset "errors: weight and value ranges, minimize" begin
        @test_throws ArgumentError solve(net, obs; soft = Dict("B" => (soft_val, 1.5)))
        @test_throws ArgumentError solve(net, obs; soft = Dict("B" => (soft_val, -0.1)))
        @test_throws ArgumentError solve(net, obs; soft = Dict("B" => (150.0, 0.5)))
        withenv("DS_SCC_METHOD" => "minimize") do
            @test_throws ArgumentError solve(net, obs; soft = Dict("B" => (soft_val, 0.5)))
        end
    end

    @testset "a cyclic component honours the blend" begin
        # B ⇄ D mutual activation fed by A: the component's soft-observed node
        # settles at its blend and the residual verdict stays honest.
        net3 = mknet(["A", "B", "D"], [edge("A", "B"), edge("B", "D"), edge("D", "B")])
        obs3 = Dict("A" => (50.0, 1.0))
        free = solve(net3, obs3)
        free_b = act(free, "B")
        target = 2.0 / 100.0
        @test free_b > target
        r = solve(net3, obs3; soft = Dict("B" => (2.0, 0.5)))
        @test r.converged
        b = act(r, "B")
        @test target < b < free_b
        # D follows the reconciled B downward
        @test act(r, "D") < act(free, "D")
    end

    @testset "FR-001: no soft observations means a bit-identical solve" begin
        plain = solve(net, obs)
        for env in ("off", "blend")
            withenv("DS_SOFT_OBS" => env) do
                r = solve(net, obs)
                @test r.node_activities == plain.node_activities
                @test r.diagnostics["soft_obs_rule"] == "off"
                @test r.diagnostics["soft_observed"] == 0
                @test isempty(r.diagnostics["soft_overridden"])
            end
        end
    end

    @testset "the LM fallback solves the blended system (research.md R1)" begin
        # The specs/066 Goodwin ring: P and nine subunits AND-assemble into a
        # chain ending at D, which depletes every subunit — gain 9 across nine
        # lags, so the fixed point repels damped relaxation and only the LM
        # fallback converges it. A soft observation inside the component must
        # go through the SAME machinery: rows become g(x) - x, the solve
        # converges, and the answer is the cap-independent blended fixed point.
        subs = ["S$i" for i in 1:9]
        names = vcat(["P"], subs, ["C1", "C2", "C3", "D"], ["rA", "rB", "rC", "rD"])
        nodes = Dict(u => DeltaSignal.NetworkNode(u, "R-$u", "protein", nothing, u, BL)
                     for u in names)
        le(s, t, et = "input") =
            DeltaSignal.LogicNetworkEdge(s, t, et != "output", et != "depletion", 1.0, et)
        edges = [le("P", "rA"); [le(s, "rA") for s in subs]; le("rA", "C1", "output");
                 le("C1", "rB"); le("rB", "C2", "output"); le("C2", "rC");
                 le("rC", "C3", "output"); le("C3", "rD"); le("rD", "D", "output");
                 [le("D", s, "depletion") for s in subs]]
        ringnet = DeltaSignal.ReactionNetwork(nodes, edges,
            Dict{String, DeltaSignal.SetExpansionMapping}(), Set{String}(),
            Dict{String, Set{String}}(), nothing, nothing, nothing, nothing, nothing)
        rsolve(obs; soft = Dict{String, Tuple{Float64, Float64}}(), iters = 500) =
            withenv("DS_DEPLETION_OWN_PRODUCT" => "full") do
                DeltaSignal.solve_steady_state(ringnet, obs,
                    DeltaSignal.SteadyStateParams(1.0, 0.1, iters, 1e-6, "penalty");
                    soft_observations = soft)
            end
        pobs = Dict("P" => (80.0, 1.0))

        free = rsolve(pobs)
        @test free.converged            # the fixture needs and gets the fallback
        @test free.diagnostics["nonconv_components"] == 1
        free_d = act(free, "D")

        soft_y = 2.0                    # UI: 2x, far from the ring's own answer
        for w in (0.25, 0.5, 0.75)
            r = rsolve(pobs; soft = Dict("D" => (soft_y, w)))
            # the whole point: a soft-carrying component CONVERGES through LM
            @test r.converged
            @test r.diagnostics["soft_observed"] == 1
            @test min(free_d, soft_y / 100.0) < act(r, "D") < max(free_d, soft_y / 100.0)
            # cap-independent: the same blended fixed point from a longer sweep
            r2 = rsolve(pobs; soft = Dict("D" => (soft_y, w)), iters = 2000)
            @test act(r2, "D") ≈ act(r, "D") rtol = 1e-5
        end
        # weight 1 through the fallback equals the hard pin on the same node
        hard = rsolve(merge(pobs, Dict("D" => (soft_y, 1.0))))
        w1 = rsolve(pobs; soft = Dict("D" => (soft_y, 1.0)))
        @test w1.converged
        @test act(w1, "D") ≈ act(hard, "D") atol = 1e-9
        @test act(w1, "S1") ≈ act(hard, "S1") rtol = 1e-6
    end

    @testset "DS_SOFT_OBS validation" begin
        # the DEFAULT is off; set, not assumed: the serving .env passes
        # DS_SOFT_OBS=blend into the dev container (review 2026-11-20)
        withenv("DS_SOFT_OBS" => nothing) do
            @test DeltaSignal.soft_obs_mode() == "off"
        end
        withenv("DS_SOFT_OBS" => "blend") do
            @test DeltaSignal.soft_obs_mode() == "blend"
        end
        withenv("DS_SOFT_OBS" => "sideways") do
            @test_throws ArgumentError DeltaSignal.soft_obs_mode()
        end
    end
end

println("test_soft_observations.jl passed")
