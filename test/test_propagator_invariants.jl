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
        # numbers describe. Without this the suite fails with bare assertion
        # errors and reads as
        # "the propagator broke" rather than "DS_AND_MODE is overridden in this
        # shell". That is the exact confusion CLAUDE.md records from the six
        # days docker-compose.dev.yml drifted from the code defaults.
        #
        # This deliberately does NOT set the environment itself: these are
        # characterization tests for the SHIPPED default, so an accidental
        # change to that default should surface here, not be masked.
        cfg = DeltaSignal.resolve_reaction_eval_config()
        @test cfg.and_mode == "hill_sat"
        @test cfg.inhibition_mode == "divide"
        @test cfg.or_mode == "mean"
        @test cfg.assembly_limiting
        if cfg.and_mode != "hill_sat" || cfg.inhibition_mode != "divide"
            @warn "Propagator invariants describe and_mode=hill_sat / " *
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

    @testset "aggregation mode IS a no-op for a single input" begin
        # A lone activator has nothing to combine with, so AND and OR must
        # agree. Under the old hill_log default they did NOT: AND applied tanh
        # compression to a single input, giving 74.07 for a 100x source where
        # OR gave 100.0, and 0.0007 for a knockout where OR gave 0. That
        # asymmetry is gone -- specs/010.
        and_net = mknet(["P", "T"], [edge("P", "T"; is_and=true)])
        or_net  = mknet(["P", "T"], [edge("P", "T"; is_and=false)])

        for f in (0.0, 0.5, 1.0, 2.0, 100.0)
            a = fold(and_net, Dict("P" => (f, 1.0)), "T")
            o = fold(or_net,  Dict("P" => (f, 1.0)), "T")
            @test a ≈ o rtol=1e-6
        end
        # And both are the identity, which is what a single input means.
        @test fold(and_net, Dict("P" => (100.0, 1.0)), "T") ≈ 100.0 rtol=1e-6
        @test fold(and_net, Dict("P" => (0.0, 1.0)), "T") == 0.0
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

    @testset "magnitude is preserved through path depth" begin
        # A chain of SINGLE-input reactions has nothing to combine at any step,
        # so it must be the identity at every depth. Under the old hill_log
        # default it was not: tanh compression compounded, and a 100x source
        # read 74x at one hop, 33x at five and 19x at ten. Two readouts with
        # identical biology and different path lengths therefore received
        # different predicted folds, which is what made those magnitudes
        # unusable as a scale. specs/010 has the table.
        for f in (0.1, 0.5, 2.0, 10.0, 100.0)
            for d in (1, 5, 10, 20)
                @test chain_fold(d, f) ≈ f rtol=1e-6
            end
        end

        # Baseline is the fixed point, at any depth.
        for d in (1, 5, 20, 40)
            @test chain_fold(d, 1.0) ≈ 1.0 rtol=1e-9
        end

        # Depth-invariance stated directly: the same source fold gives the same
        # answer one hop away and ten.
        @test chain_fold(1, 100.0) ≈ chain_fold(10, 100.0) rtol=1e-6
        @test chain_fold(1, 0.1) ≈ chain_fold(10, 0.1) rtol=1e-6
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

    @testset "a composition edge behaves exactly like an assembly edge" begin
        # LNG emits `composition` from a complex node to the complex that
        # contains it (Reactome hasComponent, <=2 hops). The solver must treat
        # it as assembly: same limiting-reactant aggregation, bit-identical.
        N(id) = DeltaSignal.NetworkNode(id, "R-HSA-" * id, "unknown", nothing, id, 0.01)
        function comp_net(ty)
            nodes = Dict(id => N(id) for id in ("A", "B", "C", "X"))
            edges = [DeltaSignal.LogicNetworkEdge("A", "X", true, true, 1.0, ty),
                     DeltaSignal.LogicNetworkEdge("B", "X", true, true, 1.0, ty),
                     DeltaSignal.LogicNetworkEdge("C", "X", true, true, 1.0, "catalyst")]
            DeltaSignal.ReactionNetwork(nodes, edges, Dict{String,DeltaSignal.SetExpansionMapping}())
        end
        P = DeltaSignal.SteadyStateParams(1.0, 1e-6, 200, 1e-8, "penalty")
        for obs in (Dict("A" => (0.0, 1.0)), Dict("A" => (50.0, 1.0)), Dict("B" => (3.0, 1.0), "C" => (0.5, 1.0)))
            a = DeltaSignal.solve_steady_state(comp_net("assembly"), obs, P).node_activities
            c = DeltaSignal.solve_steady_state(comp_net("composition"), obs, P).node_activities
            @test a == c
        end
        # and NOT like a plain input: under DS_ASSEMBLY_LIMITING the scarcest
        # component caps the container, a plain AND input does not.
        # (A=0 would make both sides 0 under hill_sat's 0 x anything = 0 and
        # the assertion vacuous; A=0.5 discriminates: 0.005 vs 0.25)
        obs = Dict("A" => (0.5, 1.0), "B" => (50.0, 1.0))
        c = DeltaSignal.solve_steady_state(comp_net("composition"), obs, P).node_activities["X"]
        i = DeltaSignal.solve_steady_state(comp_net("input"), obs, P).node_activities["X"]
        @test c < i
    end

    @testset "DS_COMPOSITION_MODE=limit: a component can lower its container, never raise it" begin
        # specs/016. Under the default ("assembly") a composition edge joins the
        # AND cluster and meets the container's producing reaction through
        # DS_OR_COMBINE: `max` lets a baseline component MASK a DOWN arriving
        # through the producing reaction, `gate` lets an over-expressed
        # component MULTIPLY into every container above it. Both measured
        # (-79 / -79 held-out). "limit" makes the edge a pure limiter.
        N(id) = DeltaSignal.NetworkNode(id, "R-HSA-" * id, "unknown", nothing, id, 0.01)
        # P --(input)--> R --(output)--> X   the container's producing route
        # C --(composition)--> X            a component complex of X
        function lim_net(; with_comp=true)
            nodes = Dict(id => N(id) for id in ("P", "R", "X", "C"))
            edges = [DeltaSignal.LogicNetworkEdge("P", "R", true, true, 1.0, "input"),
                     DeltaSignal.LogicNetworkEdge("R", "X", false, true, 1.0, "output")]
            with_comp && push!(edges, DeltaSignal.LogicNetworkEdge("C", "X", true, true, 1.0, "composition"))
            DeltaSignal.ReactionNetwork(nodes, edges, Dict{String,DeltaSignal.SetExpansionMapping}())
        end
        P = DeltaSignal.SteadyStateParams(1.0, 1e-6, 200, 1e-8, "penalty")
        act(net, obs, env) = withenv("DS_COMPOSITION_MODE" => env) do
            DeltaSignal.solve_steady_state(net, obs, P).node_activities
        end
        down_route = Dict("P" => (0.3, 1.0))                  # producing route DOWN, component untouched
        ko_comp    = Dict("C" => (0.0, 1.0))                  # component knocked out
        half_comp  = Dict("C" => (0.5, 1.0))
        oe_comp    = Dict("C" => (80.0, 1.0))                 # component over-expressed

        # default is byte-identical to the flag being absent
        @test act(lim_net(), down_route, nothing) == act(lim_net(), down_route, "assembly")

        # the masking defect, pinned: default `max` reads a DOWN route as baseline
        # when a baseline component sits beside it; limit lets the DOWN through
        @test act(lim_net(), down_route, "assembly")["X"] / BL ≈ 1.0 rtol=1e-9
        @test act(lim_net(), down_route, "limit")["X"] ≈ act(lim_net(with_comp=false), down_route, nothing)["X"] rtol=1e-9
        @test act(lim_net(), down_route, "limit")["X"] / BL < 0.85

        # a component can only lower the container ...
        @test act(lim_net(), ko_comp, "limit")["X"] < 1e-9
        @test act(lim_net(), half_comp, "limit")["X"] / BL ≈ 0.5 rtol=1e-9
        # ... never raise it: over-expression changes nothing
        @test act(lim_net(), oe_comp, "limit")["X"] == act(lim_net(with_comp=false), Dict{String,Tuple{Float64,Float64}}(), nothing)["X"]
        # and the limiter multiplies the producing route, it does not replace it
        both = Dict("P" => (0.5, 1.0), "C" => (0.5, 1.0))
        @test act(lim_net(), both, "limit")["X"] / BL ≈ 0.25 rtol=1e-6

        # a container with NO producing route (the severed Interferon branch)
        # reads baseline x limiter, so a component knockout still reaches it
        nodes = Dict(id => N(id) for id in ("C", "X"))
        sev = DeltaSignal.ReactionNetwork(nodes, [DeltaSignal.LogicNetworkEdge("C", "X", true, true, 1.0, "composition")],
                                          Dict{String,DeltaSignal.SetExpansionMapping}())
        @test act(sev, half_comp, "limit")["X"] / BL ≈ 0.5 rtol=1e-9
        @test act(sev, oe_comp, "limit")["X"] / BL ≈ 1.0 rtol=1e-9

        # edge order does not matter
        rev = let n = lim_net(); DeltaSignal.ReactionNetwork(n.nodes, reverse(n.edges), n.set_mappings) end
        @test act(rev, both, "limit") == act(lim_net(), both, "limit")

        # a typo is an error, not a silent model choice
        @test_throws ArgumentError act(lim_net(), both, "Limit")
    end

    @testset "DS_COMPOSITION_MODE=limit_novel: a hierarchy edge that repeats a producing route is skipped" begin
        # specs/016. 54% of composition edges connect a component to a container
        # that a reaction ALREADY builds from it. Under `limit` the component's
        # fold enters twice (reaction input and hierarchy edge) and squares at
        # every level of DSB Repair's 35-deep nesting -- RAD52 OE zeroed 281 of
        # 365 composition sources. `limit_novel` keeps only the hops with no
        # reaction, which are the ones that carried Interferon alpha/beta.
        N(id) = DeltaSignal.NetworkNode(id, "R-HSA-" * id, "unknown", nothing, id, 0.01)
        # C --(input)--> R --(output)--> X   and   C --(composition)--> X   (redundant)
        # C --(composition)--> Y            with no reaction building Y   (novel)
        function nets(; comp=true)
            nodes = Dict(id => N(id) for id in ("C", "R", "X", "Y"))
            edges = [DeltaSignal.LogicNetworkEdge("C", "R", true, true, 1.0, "input"),
                     DeltaSignal.LogicNetworkEdge("R", "X", false, true, 1.0, "output")]
            if comp
                push!(edges, DeltaSignal.LogicNetworkEdge("C", "X", true, true, 1.0, "composition"))
                push!(edges, DeltaSignal.LogicNetworkEdge("C", "Y", true, true, 1.0, "composition"))
            end
            DeltaSignal.ReactionNetwork(nodes, edges, Dict{String,DeltaSignal.SetExpansionMapping}())
        end
        P = DeltaSignal.SteadyStateParams(1.0, 1e-6, 200, 1e-8, "penalty")
        act(net, obs, env) = withenv("DS_COMPOSITION_MODE" => env) do
            DeltaSignal.solve_steady_state(net, obs, P).node_activities
        end
        half = Dict("C" => (0.5, 1.0))
        # under `limit` the redundant edge squares the fold at X
        @test act(nets(), half, "limit")["X"] / BL ≈ 0.25 rtol=1e-6
        # under `limit_novel` X reads the producing route only: identical to no composition edge
        @test act(nets(), half, "limit_novel")["X"] ≈ act(nets(comp=false), half, nothing)["X"] rtol=1e-9
        @test act(nets(), half, "limit_novel")["X"] / BL ≈ 0.5 rtol=1e-6
        # and the novel hop still limits: Y has no producing reaction and follows its component down
        @test act(nets(), half, "limit_novel")["Y"] / BL ≈ 0.5 rtol=1e-6
        @test act(nets(), Dict("C" => (80.0, 1.0)), "limit_novel")["Y"] / BL ≈ 1.0 rtol=1e-9
        # default untouched, edge order irrelevant
        @test act(nets(), half, nothing) == act(nets(), half, "assembly")
        rev = let n = nets(); DeltaSignal.ReactionNetwork(n.nodes, reverse(n.edges), n.set_mappings) end
        @test act(rev, half, "limit_novel") == act(nets(), half, "limit_novel")
    end

    @testset "DS_DEPLETION_OWN_PRODUCT: a product that is low because its substrate is low restores nothing" begin
        # specs/016, traced on AKT1-KO -> TP53: nuclear p-MDM2 is depleted by
        # the dimers and complexes built FROM it. When its supply halves they
        # halve too, and the depletion edge de-represses it back (0.50 -> 0.82).
        N(id) = DeltaSignal.NetworkNode(id, "R-HSA-" * id, "unknown", nothing, id, 0.01)
        # U -> R0 -> X (supply);  X + W -> R -> P (P made from X and partner W);  P -| X (depletion)
        function dep_net(; depleter_is_product=true)
            ids = ("U", "R0", "X", "W", "R", "P", "Q")
            nodes = Dict(id => N(id) for id in ids)
            edges = [DeltaSignal.LogicNetworkEdge("U", "R0", true, true, 1.0, "input"),
                     DeltaSignal.LogicNetworkEdge("R0", "X", false, true, 1.0, "output"),
                     DeltaSignal.LogicNetworkEdge("X", "R", true, true, 1.0, "input"),
                     DeltaSignal.LogicNetworkEdge("W", "R", true, true, 1.0, "input"),
                     DeltaSignal.LogicNetworkEdge("R", "P", false, true, 1.0, "output"),
                     DeltaSignal.LogicNetworkEdge(depleter_is_product ? "P" : "Q", "X", true, false, 1.0, "depletion")]
            DeltaSignal.ReactionNetwork(nodes, edges, Dict{String,DeltaSignal.SetExpansionMapping}())
        end
        P = DeltaSignal.SteadyStateParams(1.0, 1e-6, 500, 1e-8, "penalty")
        act(net, obs, env) = withenv("DS_DEPLETION_OWN_PRODUCT" => env) do
            DeltaSignal.solve_steady_state(net, obs, P).node_activities
        end
        half_supply = Dict("U" => (0.5, 1.0))
        partner_ko  = Dict("W" => (0.0, 1.0))
        partner_oe  = Dict("W" => (20.0, 1.0))

        # default is byte-identical to the flag being absent
        @test act(dep_net(), half_supply, nothing) == act(dep_net(), half_supply, "full")

        # the defect, pinned: under "full" a halved supply reads back above 0.5 at X
        @test act(dep_net(), half_supply, "full")["X"] / BL > 0.55
        # suppress_only: the supply drop passes through undiluted
        @test act(dep_net(), half_supply, "suppress_only")["X"] / BL ≈ 0.5 rtol=1e-6
        @test act(dep_net(), half_supply, "suppress_only")["P"] / BL < act(dep_net(), half_supply, "full")["P"] / BL

        # suppression by an abundant own product is kept: partner OE lifts P and drains X in both modes
        @test act(dep_net(), partner_oe, "suppress_only")["X"] / BL < 0.85
        @test act(dep_net(), partner_oe, "suppress_only")["X"] == act(dep_net(), partner_oe, "full")["X"]

        # what is deliberately given up: partner KO no longer de-represses X through P
        @test act(dep_net(), partner_ko, "full")["X"] / BL > 1.15
        @test act(dep_net(), partner_ko, "suppress_only")["X"] / BL ≈ 1.0 rtol=1e-9

        # a depleter that is NOT the target's product is untouched by the mode
        foreign = dep_net(depleter_is_product=false)
        for obs in (Dict("Q" => (0.0, 1.0)), Dict("Q" => (5.0, 1.0)))
            @test act(foreign, obs, "suppress_only") == act(foreign, obs, "full")
        end
        @test act(foreign, Dict("Q" => (0.0, 1.0)), "suppress_only")["X"] / BL > 1.15

        # edge order does not matter, and a typo is an error
        rev = let n = dep_net(); DeltaSignal.ReactionNetwork(n.nodes, reverse(n.edges), n.set_mappings) end
        @test act(rev, half_supply, "suppress_only") == act(dep_net(), half_supply, "suppress_only")
        @test_throws ArgumentError act(dep_net(), half_supply, "suppress-only")
    end

    @testset "review fixes 2026-09-19: dedup keeps the input role, own-product ignores catalysts, elasticity acts under the flat solve" begin
        N(id) = DeltaSignal.NetworkNode(id, "R-HSA-" * id, "unknown", nothing, id, 0.01)
        P = DeltaSignal.SteadyStateParams(1.0, 1e-6, 300, 1e-8, "penalty")
        mk(nodes, edges) = DeltaSignal.ReactionNetwork(Dict(id => N(id) for id in nodes), edges, Dict{String,DeltaSignal.SetExpansionMapping}())
        E(s, t, is_and, pos, ty) = DeltaSignal.LogicNetworkEdge(s, t, is_and, pos, 1.0, ty)

        # (a) S feeds X as a plain input AND by a composition edge. Under dedup the
        # plain role must win: a 4x input stays a 4x input, not a <=1 limiter.
        dup = mk(("S", "X", "K"), [E("S", "X", true, true, "input"), E("S", "X", true, true, "composition"), E("K", "X", true, true, "catalyst")])
        x_on  = withenv("DS_DEDUP_ACTIVATORS" => "1", "DS_COMPOSITION_MODE" => "limit") do
            DeltaSignal.solve_steady_state(dup, Dict("S" => (4.0, 1.0)), P).node_activities["X"] end
        x_off = withenv("DS_DEDUP_ACTIVATORS" => "0", "DS_COMPOSITION_MODE" => "limit") do
            DeltaSignal.solve_steady_state(dup, Dict("S" => (4.0, 1.0)), P).node_activities["X"] end
        @test x_on / BL > 1.15
        @test x_on ≈ x_off rtol=1e-9

        # (b) X is the CATALYST of R (S + X -> R -> P) and P depletes X. P is not
        # X's own product; partner KO must still de-repress X under suppress_only.
        cat = mk(("S", "X", "R", "P"), [E("S", "R", true, true, "input"), E("X", "R", true, true, "catalyst"),
                                        E("R", "P", false, true, "output"), E("P", "X", true, false, "depletion")])
        for mode in ("full", "suppress_only")
            x = withenv("DS_DEPLETION_OWN_PRODUCT" => mode) do
                DeltaSignal.solve_steady_state(cat, Dict("S" => (0.0, 1.0)), P).node_activities["X"] end
            @test x / BL > 1.15
        end

        # (d) elasticity needs SCC membership; with the legacy flat solve it used to
        # be a silent no-op. A positive two-node loop driven above baseline must
        # now respond to eps.
        loop = mk(("U", "A", "B"), [E("U", "A", true, true, "input"), E("A", "B", true, true, "input"), E("B", "A", true, true, "input")])
        a1 = withenv("DS_SCC_SOLVE" => "0", "DS_LOOP_ELASTICITY" => nothing) do
            DeltaSignal.solve_steady_state(loop, Dict("U" => (1.1, 1.0)), P).node_activities["A"] end
        a5 = withenv("DS_SCC_SOLVE" => "0", "DS_LOOP_ELASTICITY" => "0.5") do
            DeltaSignal.solve_steady_state(loop, Dict("U" => (1.1, 1.0)), P).node_activities["A"] end
        @test a1 != a5
        @test a5 / BL < a1 / BL
    end

    @testset "DS_COMPOSITION_GROUP: shards of one entity are alternatives, components are co-required" begin
        # HDR case: BCDX2 exists as 33 node copies (one per variant reaction) and
        # every copy feeds the container by a composition edge. Under plain
        # assembly-limiting the container is capped by the LEAST copy, so a
        # perturbation reaching one variant reaction collapses it while 32 copies
        # sit at baseline. Grouping composition inputs by base entity -- max
        # within a group, min across groups -- makes copies of one entity
        # alternatives and distinct components co-required.
        N(id, stid) = DeltaSignal.NetworkNode(id, stid, "unknown", nothing, id, 0.01)
        # two shards of BCDX2 (same stable id, one a ::variant), one CX3, one container
        nodes = Dict("b1" => N("b1", "R-HSA-5685316"), "b2" => N("b2", "R-HSA-5685316::variant::R-HSA-1"),
                     "c1" => N("c1", "R-HSA-5685311"), "X" => N("X", "R-HSA-5686104"))
        edges = [DeltaSignal.LogicNetworkEdge(s, "X", true, true, 1.0, "composition") for s in ("b1", "b2", "c1")]
        net = DeltaSignal.ReactionNetwork(nodes, edges, Dict{String,DeltaSignal.SetExpansionMapping}())
        P = DeltaSignal.SteadyStateParams(1.0, 1e-6, 200, 1e-8, "penalty")
        ko_one_shard = Dict("b1" => (0.0, 1.0))                  # one BCDX2 copy knocked out, the other at baseline
        ko_component = Dict("c1" => (0.0, 1.0))                  # CX3 knocked out
        up_one_shard = Dict("b1" => (50.0, 1.0))

        # default (grouping off): byte-identical to the flag being absent
        ref = withenv("DS_COMPOSITION_GROUP" => nothing) do
            DeltaSignal.solve_steady_state(net, ko_one_shard, P)
        end
        off = withenv("DS_COMPOSITION_GROUP" => "0") do
            DeltaSignal.solve_steady_state(net, ko_one_shard, P)
        end
        @test off.node_activities == ref.node_activities
        # and under the default the least copy caps the container: one shard KO collapses X
        @test off.node_activities["X"] < 1e-6

        # grouping ON: the surviving copy carries the entity, X stays at baseline
        on = withenv("DS_COMPOSITION_GROUP" => "1") do
            DeltaSignal.solve_steady_state(net, ko_one_shard, P)
        end
        @test isapprox(on.node_activities["X"], 0.01; rtol = 1e-6)
        # a distinct component knocked out still collapses X (co-required)
        onc = withenv("DS_COMPOSITION_GROUP" => "1") do
            DeltaSignal.solve_steady_state(net, ko_component, P)
        end
        @test onc.node_activities["X"] < 1e-6
        # one copy UP with grouping: max within the group lifts X (min across groups then caps at CX3 = baseline)
        onu = withenv("DS_COMPOSITION_GROUP" => "1") do
            DeltaSignal.solve_steady_state(net, up_one_shard, P)
        end
        @test isapprox(onu.node_activities["X"], 0.01; rtol = 1e-6)   # capped by CX3 at baseline
        # a typo throws rather than silently choosing (spellings no bool parser accepts)
        for bad in ("maybe", "1.5", "tru")
            withenv("DS_COMPOSITION_GROUP" => bad) do
                @test_throws ArgumentError DeltaSignal.solve_steady_state(net, ko_one_shard, P)
            end
        end
    end
end

