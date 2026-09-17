#!/usr/bin/env julia
#
# Propagator invariants, and the asymmetries between the activation and
# inhibition paths.
#
# These are CHARACTERIZATION tests: they pin what the propagator does today so
# a change to it is visible, not a claim that today's behaviour is right. Two
# of the facts below are surprising enough to be worth stating out loud.
#
# Found while sweeping solver invariants for the false-change error mode. The
# hypotheses that DIED are as useful as the ones that held, so they are
# asserted here too rather than described in a commit message nobody reads:
# monotonicity holds in both directions, and nothing drifts off baseline.
#
# NOTE ON FIXTURES: LogicNetworkEdge is (parent, child, is_and, is_positive,
# stoichiometry, edge_type). is_and comes BEFORE is_positive. Getting them the
# wrong way round silently turns an OR activator into an inhibitor, and the
# resulting numbers look like an inverted model rather than a broken fixture.

using Test

include(joinpath(@__DIR__, "..", "src", "DeltaSignal.jl"))
using .DeltaSignal
using .DeltaSignal: resolve_reaction_eval_config

const BL = 0.01
const PARAMS = DeltaSignal.SteadyStateParams(1.0, 0.1, 500, 1e-6, "penalty")

mknet(nodes, edges) = DeltaSignal.ReactionNetwork(
    Dict(n => DeltaSignal.NetworkNode(n, n, "protein", nothing, n, BL) for n in nodes),
    edges, Dict{String, DeltaSignal.SetExpansionMapping}())

edge(s, t; is_and::Bool=true, positive::Bool=true) =
    DeltaSignal.LogicNetworkEdge(s, t, is_and, positive, 1.0, "input")

function fold(net, obs, node)
    r = DeltaSignal.solve_steady_state(net, obs, PARAMS)
    return r.node_activities[node] / BL
end

"""A chain of `depth` single-input reactions, source pinned at `f`."""
function chain_fold(depth::Int, f::Float64)
    nodes = ["N$i" for i in 0:depth]
    net = mknet(nodes, [edge("N$(i-1)", "N$i") for i in 1:depth])
    return fold(net, Dict("N0" => (f, 1.0)), "N$depth")
end

@testset "propagator invariants" begin

    @testset "the config these numbers describe" begin
        # Fail FAST and by name if the resolved config is not the one these
        # magic numbers were measured under. Without this the suite fails with
        # bare assertion errors -- 74.07 became something else -- and reads as
        # "the propagator broke" rather than "DS_AND_MODE is overridden in this
        # shell". That is the exact confusion CLAUDE.md records from the six
        # days docker-compose.dev.yml drifted from the code defaults.
        #
        # This deliberately does NOT set the environment itself: these are
        # characterization tests for the SHIPPED default, so an accidental
        # change to that default should surface here, not be masked.
        cfg = DeltaSignal.resolve_reaction_eval_config()
        @test cfg.and_mode == "hill_log"
        @test cfg.inhibition_mode == "divide"
        @test cfg.or_mode == "mean"
        @test cfg.assembly_limiting
        if cfg.and_mode != "hill_log" || cfg.inhibition_mode != "divide"
            @warn "Propagator invariants describe and_mode=hill_log / " *
                  "inhibition_mode=divide. The resolved config differs, so the " *
                  "magnitudes below WILL fail. Check DS_* in this environment." cfg
        end
    end

    @testset "inhibition is exactly reciprocal" begin
        # The `divide` form: H = (b+eps)/(x+eps). With one activator at
        # baseline, the target reads 1/inhibitor exactly. No compression is
        # applied on this path -- unlike the activation path below.
        net = mknet(["A", "I", "T"], [edge("A", "T"), edge("I", "T"; positive=false)])
        g(v) = fold(net, Dict("A" => (1.0, 1.0), "I" => (v, 1.0)), "T")
        @test g(0.5) ≈ 2.0    rtol=1e-6
        @test g(2.0) ≈ 0.5    rtol=1e-6
        @test g(10.0) ≈ 0.1   rtol=1e-6
        @test g(100.0) ≈ 0.01 rtol=1e-6
        @test g(1.0) ≈ 1.0    rtol=1e-9
    end

    @testset "de-repression saturates at h_max while activation does not" begin
        # Removing an inhibitor caps at DS_HILL_SAT_H_MAX (10.0). Driving an
        # activator to UI 100 reaches 74.07. The two directions therefore have
        # different ceilings, which matters for any fold-magnitude claim.
        inh = mknet(["A", "I", "T"], [edge("A", "T"), edge("I", "T"; positive=false)])
        @test fold(inh, Dict("A" => (1.0, 1.0), "I" => (0.0, 1.0)), "T") ≈ 10.0 rtol=1e-6

        act = mknet(["P", "T"], [edge("P", "T")])
        @test fold(act, Dict("P" => (100.0, 1.0)), "T") > 10.0
    end

    @testset "aggregation mode is NOT a no-op for a single input" begin
        # A lone activator has nothing to combine with, so AND and OR should
        # agree. They do not: the AND path applies hill_log compression to a
        # single input, the OR path (mean of one) is exact.
        and_net = mknet(["P", "T"], [edge("P", "T"; is_and=true)])
        or_net  = mknet(["P", "T"], [edge("P", "T"; is_and=false)])

        @test fold(or_net,  Dict("P" => (100.0, 1.0)), "T") ≈ 100.0 rtol=1e-6
        @test fold(and_net, Dict("P" => (100.0, 1.0)), "T") ≈ 74.0673 rtol=1e-4
        @test fold(and_net, Dict("P" => (100.0, 1.0)), "T") <
              fold(or_net,  Dict("P" => (100.0, 1.0)), "T")

        # A knockout is exact through OR and very nearly exact through AND.
        @test fold(or_net,  Dict("P" => (0.0, 1.0)), "T") ≈ 0.0 atol=1e-9
        @test fold(and_net, Dict("P" => (0.0, 1.0)), "T") < 1e-3
    end

    @testset "monotonicity holds in both directions" begin
        # A non-monotonic response would make the predicted direction erratic.
        # Swept across the full operating range, including both rails.
        act = mknet(["P1", "P2", "T"], [edge("P1", "T"), edge("P2", "T")])
        prev = -Inf
        for v in (0.0, 0.1, 0.25, 0.5, 0.8, 1.0, 1.5, 2.0, 5.0, 10.0, 50.0, 100.0)
            out = fold(act, Dict("P1" => (v, 1.0), "P2" => (1.0, 1.0)), "T")
            @test out >= prev - 1e-12
            prev = out
        end

        inh = mknet(["A", "I", "T"], [edge("A", "T"), edge("I", "T"; positive=false)])
        prev = Inf
        for v in (0.0, 0.1, 0.5, 1.0, 2.0, 10.0, 100.0)
            out = fold(inh, Dict("A" => (1.0, 1.0), "I" => (v, 1.0)), "T")
            @test out <= prev + 1e-12
            prev = out
        end
    end

    @testset "signal decays toward baseline with path depth" begin
        # A single-input chain is not identity: hill_log compression compounds.
        # UI 100 arrives as 19.1 after ten hops, which is roughly the median
        # path length in the catalog.
        @test chain_fold(1, 100.0) ≈ 74.0673 rtol=1e-4
        @test chain_fold(10, 100.0) ≈ 19.14 rtol=1e-2
        @test chain_fold(40, 100.0) < chain_fold(10, 100.0)

        # Baseline is the fixed point: it must not decay at any depth.
        for d in (1, 5, 20, 40)
            @test chain_fold(d, 1.0) ≈ 1.0 rtol=1e-12
        end

        # The decay is SYMMETRIC in log-fold -- it is a function of magnitude,
        # not of direction. This is why it does not bias the model up or down,
        # and why it is not the explanation for the de-repression deficit.
        retained(f) = log2(chain_fold(10, f)) / log2(f)
        @test retained(0.01) ≈ retained(100.0) rtol=0.02
        @test retained(0.1) ≈ retained(10.0)  rtol=0.02
        @test 0.6 < retained(100.0) < 0.7
    end

    @testset "nothing drifts off baseline" begin
        # Hypotheses that died. A drift here would read as a false change,
        # which is the largest single error category, so they stay asserted.
        for n in (1, 2, 5, 20)
            inh_net = mknet(vcat(["I$i" for i in 1:n], ["T"]),
                            [edge("I$i", "T"; positive=false) for i in 1:n])
            @test fold(inh_net, Dict("I$i" => (1.0, 1.0) for i in 1:n), "T") ≈ 1.0 rtol=1e-12
        end

        # An unperturbed negative-feedback loop sits at baseline and converges.
        loop = mknet(["A", "B"], [edge("A", "B"), edge("B", "A"; positive=false)])
        r = DeltaSignal.solve_steady_state(loop, Dict{String, Tuple{Float64, Float64}}(), PARAMS)
        @test r.node_activities["A"] / BL ≈ 1.0 rtol=1e-12
        @test r.node_activities["B"] / BL ≈ 1.0 rtol=1e-12
        @test r.converged
    end
end
