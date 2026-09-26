# Loop elasticity — DS_LOOP_ELASTICITY
#
# AND is multiplication of fold-changes, so every activator edge has elasticity
# d log(out)/d log(in) = 1, and a positive feedback cycle has loop gain exactly
# 1 at baseline. Baseline is then a knife-edge, not a stable state: on the
# two-reaction AND-loop below, U = 1.01 drives A to 1.09 / 2.30 / 100.0 at
# 50 / 500 / 5000 sweeps, and U = 0.90 collapses it to 0.0026 and reports
# converged. The answer is wherever the iteration stopped, which is why sweep
# order (UUID hash) and iteration budget both change catalog predictions.
#
# With elasticity eps < 1 applied to the CYCLE-CLOSING edges only (source and
# target in one SCC), the loop equation becomes A = U * A^(eps^4) for this
# fixture (four in-loop edges), so A = U^(1/(1 - eps^4)): a unique, stable
# positive root. Zero becomes unstable (A^eps >> A near zero), so collapse
# needs a real in-loop knockout. Acyclic edges are untouched.
#
# Baseline is 0.01 internal / 1.0 UI throughout.

using Test
using DeltaSignal
const DS = DeltaSignal

N(id) = DS.NetworkNode(id, "R-HSA-" * id, "unknown", nothing, id, 0.01)
E(s, t, and, pos, ty) = DS.LogicNetworkEdge(s, t, and, pos, 1.0, ty)

"""Positive feedback, AND-closed: U -> A ; A -> r1 -> B ; B -> r2 -> A."""
function posloop()
    nodes = Dict(id => N(id) for id in ("U", "A", "B", "r1", "r2"))
    edges = [E("U", "A", true, true, "input"),
             E("A", "r1", true, true, "input"), E("r1", "B", true, true, "output"),
             E("B", "r2", true, true, "input"), E("r2", "A", true, true, "output")]
    DS.ReactionNetwork(nodes, edges, Dict{String,DS.SetExpansionMapping}())
end

"""Straight chain, no cycle: S -> r1 -> M -> r2 -> T."""
function chain()
    nodes = Dict(id => N(id) for id in ("S", "r1", "M", "r2", "T"))
    edges = [E("S", "r1", true, true, "input"), E("r1", "M", true, true, "output"),
             E("M", "r2", true, true, "input"), E("r2", "T", true, true, "output")]
    DS.ReactionNetwork(nodes, edges, Dict{String,DS.SetExpansionMapping}())
end

P(it) = DS.SteadyStateParams(1.0, 1e-6, it, 1e-6, "penalty")

solveA(net, U, it, eps) = withenv("DS_LOOP_ELASTICITY" => eps) do
    DS.solve_steady_state(net, Dict("U" => (U, 1.0)), P(it))
end

@testset "loop elasticity" begin

    @testset "DS_LOOP_ELASTICITY guard rails" begin
        net = chain()
        obs = Dict("S" => (50.0, 1.0))
        for bad in ("0", "0.0", "-0.5", "1.5", "2", "abc", "", "inf", "nan")
            withenv("DS_LOOP_ELASTICITY" => bad) do
                @test_throws ArgumentError DS.solve_steady_state(net, obs, P(100))
            end
        end
        for good in ("1.0", "1", "0.5", "0.01")
            withenv("DS_LOOP_ELASTICITY" => good) do
                @test DS.solve_steady_state(net, obs, P(100)) !== nothing
            end
        end
    end

    @testset "default (1.0) is byte-identical to unset" begin
        for net in (posloop(), chain()), U in (1.01, 0.9, 80.0)
            obsname = net === nothing ? "" : (haskey(net.nodes, "U") ? "U" : "S")
            obs = Dict(obsname => (U, 1.0))
            ref = withenv("DS_LOOP_ELASTICITY" => nothing) do
                DS.solve_steady_state(net, obs, P(500))
            end
            got = withenv("DS_LOOP_ELASTICITY" => "1.0") do
                DS.solve_steady_state(net, obs, P(500))
            end
            @test got.node_activities == ref.node_activities
        end
    end

    @testset "acyclic networks are unaffected by any elasticity" begin
        net = chain()
        for U in (0.0, 0.5, 1.0, 1.01, 50.0, 100.0), eps in ("0.5", "0.2")
            obs = Dict("S" => (U, 1.0))
            a = withenv("DS_LOOP_ELASTICITY" => "1.0") do
                DS.solve_steady_state(net, obs, P(200)) end
            b = withenv("DS_LOOP_ELASTICITY" => eps) do
                DS.solve_steady_state(net, obs, P(200)) end
            @test a.node_activities == b.node_activities
        end
    end

    @testset "the knife-edge at eps = 1 (pinned so a fix announces itself)" begin
        net = posloop()
        # Exactly at baseline the loop holds; that is the only U where it does.
        @test isapprox(solveA(net, 1.0, 5000, "1.0").node_activities["A"] * 100, 1.0; atol=1e-9)
        # A 1% nudge rails the loop given enough sweeps, and the value depends
        # on the budget. Both are the defect.
        a50   = solveA(net, 1.01, 50,   "1.0").node_activities["A"] * 100
        a500  = solveA(net, 1.01, 500,  "1.0").node_activities["A"] * 100
        a5000 = solveA(net, 1.01, 5000, "1.0").node_activities["A"] * 100
        @test a50 < a500 < a5000
        @test a5000 > 99.0
        @test_broken isapprox(a5000, 1.01; rtol=0.05)   # what a stable loop would give
        # A 10% drop collapses it and calls that converged.
        r = solveA(net, 0.9, 500, "1.0")
        @test r.node_activities["A"] * 100 < 0.01
        @test r.converged
    end

    @testset "eps < 1 makes baseline a stable root" begin
        net = posloop()
        eps = 0.5
        # Four in-loop edges => A = U^(1/(1 - eps^4)).
        expo = 1 / (1 - eps^4)
        for U in (1.01, 0.99, 1.1, 0.9, 1.5, 0.5)
            vals = [solveA(net, U, it, string(eps)).node_activities["A"] * 100
                    for it in (50, 500, 5000)]
            # Budget-independent: the answer is a property of the network now.
            @test isapprox(vals[1], vals[2]; rtol=1e-6)
            @test isapprox(vals[2], vals[3]; rtol=1e-6)
            # And it is the analytic root.
            @test isapprox(vals[3], U^expo; rtol=1e-3)
            @test solveA(net, U, 500, string(eps)).converged
        end
    end

    @testset "eps < 1 still lets a strong drive rail and a knockout collapse" begin
        net = posloop()
        for eps in ("0.5", "0.8")
            @test solveA(net, 80.0, 500, eps).node_activities["A"] * 100 > 99.0
            @test solveA(net, 0.0, 500, eps).node_activities["A"] * 100 < 1e-6
        end
    end

    @testset "elasticity is monotone: closer to 1 means closer to the knife-edge" begin
        net = posloop()
        # For U > 1 the stable root rises toward the rail as eps -> 1.
        roots = [solveA(net, 1.1, 5000, string(e)).node_activities["A"] * 100
                 for e in (0.3, 0.5, 0.7, 0.9)]
        @test issorted(roots)
        @test roots[1] < roots[end]
    end

    # ---------------------------------------------------------------------
    # Sigmoidal elasticity (DS_LOOP_ELASTICITY_WIDTH). Constant eps < 1 was
    # measured to trade false change for MISSED change ~1:1 on the catalog:
    # 98% of the cases it broke were genuine changes flattened to NORMAL,
    # median 1.56x. In a giant SCC nearly every edge is "in-loop", so a
    # constant eps compresses every path. The leaks to damp sit at ~1.0x; the
    # signals to keep sit at >= 1.5x -- separable by magnitude.
    #
    #     eps(f) = eps_lo + (1 - eps_lo) * tanh(|log f| / w)
    #
    # Inside the band (|log f| << w) the loop is stable; outside it the edge
    # reads fold^1 again and a positive loop becomes a genuine SWITCH: an
    # above-threshold signal rails, a below-threshold leak decays.
    # ---------------------------------------------------------------------

    solveAw(net, U, it, eps, w) = withenv("DS_LOOP_ELASTICITY" => eps,
                                        "DS_LOOP_ELASTICITY_WIDTH" => w) do
        DS.solve_steady_state(net, Dict("U" => (U, 1.0)), P(it))
    end

    @testset "DS_LOOP_ELASTICITY_WIDTH guard rails" begin
        net = chain(); obs = Dict("S" => (50.0, 1.0))
        for bad in ("-1", "-0.1", "abc", "inf", "nan")
            withenv("DS_LOOP_ELASTICITY" => "0.5", "DS_LOOP_ELASTICITY_WIDTH" => bad) do
                @test_throws ArgumentError DS.solve_steady_state(net, obs, P(100))
            end
        end
        for good in ("0", "0.0", "0.3", "1", "10")
            withenv("DS_LOOP_ELASTICITY" => "0.5", "DS_LOOP_ELASTICITY_WIDTH" => good) do
                @test DS.solve_steady_state(net, obs, P(100)) !== nothing
            end
        end
    end

    @testset "width 0 is byte-identical to constant elasticity" begin
        net = posloop()
        for U in (1.01, 0.9, 1.5, 3.0, 0.2), eps in ("0.5", "0.8")
            a = solveA(net, U, 500, eps).node_activities
            b = solveAw(net, U, 500, eps, "0").node_activities
            @test a == b
        end
    end

    @testset "sigmoid: leaks inside the band decay, at every width" begin
        net = posloop()
        for w in ("0.1", "0.3", "1.0"), U in (1.01, 0.99)
            a = solveAw(net, U, 5000, "0.5", w).node_activities["A"] * 100
            @test isapprox(a, U; rtol=2e-3)            # held near U, not railed
            # and budget-independent
            @test isapprox(a, solveAw(net, U, 50, "0.5", w).node_activities["A"] * 100; rtol=1e-6)
        end
    end

    @testset "sigmoid: the loop is a switch with threshold set by the width" begin
        net = posloop()
        # w = 0.3: a 10% signal is inside the band and held; 1.5x is outside and rails.
        @test isapprox(solveAw(net, 1.1, 2000, "0.5", "0.3").node_activities["A"] * 100, 1.1; rtol=0.05)
        @test solveAw(net, 1.5, 2000, "0.5", "0.3").node_activities["A"] * 100 > 99.0
        # w = 1.0: 1.5x is inside and held; 2x is outside and rails.
        @test isapprox(solveAw(net, 1.5, 2000, "0.5", "1.0").node_activities["A"] * 100, 1.5; rtol=0.1)
        @test solveAw(net, 2.0, 2000, "0.5", "1.0").node_activities["A"] * 100 > 99.0
        # Downward: a 2x drop is outside every band tested and collapses.
        for w in ("0.3", "1.0")
            @test solveAw(net, 0.5, 2000, "0.5", w).node_activities["A"] * 100 < 1e-2
        end
        # Wider band => higher threshold, monotone.
        r03 = solveAw(net, 1.5, 2000, "0.5", "0.3").node_activities["A"]
        r10 = solveAw(net, 1.5, 2000, "0.5", "1.0").node_activities["A"]
        @test r03 > r10
    end

    @testset "sigmoid leaves acyclic networks untouched" begin
        net = chain()
        for U in (0.0, 0.9, 1.01, 1.5, 50.0), w in ("0.3", "1.0")
            a = withenv("DS_LOOP_ELASTICITY" => "1.0") do
                DS.solve_steady_state(net, Dict("S" => (U, 1.0)), P(200)) end
            b = withenv("DS_LOOP_ELASTICITY" => "0.5", "DS_LOOP_ELASTICITY_WIDTH" => w) do
                DS.solve_steady_state(net, Dict("S" => (U, 1.0)), P(200)) end
            @test a.node_activities == b.node_activities
        end
    end

    # ---------------------------------------------------------------------
    # Sigmoid ceiling (DS_LOOP_ELASTICITY_HI). With the ceiling at 1.0 an
    # above-band signal is read at fold^1, restoring loop gain 1 -- and the
    # knife-edge -- for every signal outside the band. Measured on the catalog:
    # relabelling UUIDs then moved 73 predictions (vs 14 for the original
    # solver and 0 for constant eps). With eps_hi < 1 the loop's gain stays
    # below 1 at every amplitude, so the root is unique everywhere: on this
    # fixture an above-band signal settles at U^(1/(1 - eps_hi^4)) instead of
    # railing. Direction is preserved, which is all a 3-class score needs.
    # ---------------------------------------------------------------------

    solveAh(net, U, it, hi) = withenv("DS_LOOP_ELASTICITY" => "0.5",
                                     "DS_LOOP_ELASTICITY_WIDTH" => "0.2",
                                     "DS_LOOP_ELASTICITY_HI" => hi) do
        DS.solve_steady_state(net, Dict("U" => (U, 1.0)), P(it))
    end

    @testset "DS_LOOP_ELASTICITY_HI guard rails" begin
        net = chain(); obs = Dict("S" => (50.0, 1.0))
        for bad in ("0", "-1", "1.5", "abc", "inf", "nan")
            withenv("DS_LOOP_ELASTICITY" => "0.5", "DS_LOOP_ELASTICITY_WIDTH" => "0.2",
                    "DS_LOOP_ELASTICITY_HI" => bad) do
                @test_throws ArgumentError DS.solve_steady_state(net, obs, P(100))
            end
        end
        for good in ("1", "1.0", "0.95", "0.9", "0.5")
            withenv("DS_LOOP_ELASTICITY" => "0.5", "DS_LOOP_ELASTICITY_WIDTH" => "0.2",
                    "DS_LOOP_ELASTICITY_HI" => good) do
                @test DS.solve_steady_state(net, obs, P(100)) !== nothing
            end
        end
    end

    @testset "ceiling 1.0 is byte-identical to unset" begin
        net = posloop()
        for U in (1.01, 1.5, 3.0, 0.5)
            ref = withenv("DS_LOOP_ELASTICITY" => "0.5", "DS_LOOP_ELASTICITY_WIDTH" => "0.2",
                          "DS_LOOP_ELASTICITY_HI" => nothing) do
                DS.solve_steady_state(net, Dict("U" => (U, 1.0)), P(500)) end
            @test solveAh(net, U, 500, "1.0").node_activities == ref.node_activities
        end
    end

    @testset "ceiling < 1: an above-band signal settles at a finite root" begin
        net = posloop()
        for hi in (0.9, 0.95), U in (1.5, 2.0)
            expo = 1 / (1 - hi^4)
            a = solveAh(net, U, 5000, string(hi)).node_activities["A"] * 100
            if U^expo < 99.0
                @test isapprox(a, U^expo; rtol=2e-2)     # the analytic root, not the rail
                @test a > 1.15                             # and still classified UP
                # Converged to ONE root: 500 and 5000 sweeps agree.
                @test isapprox(a, solveAh(net, U, 500, string(hi)).node_activities["A"] * 100; rtol=1e-3)
            else
                @test a > 99.0                             # 2.0^5.4 = 42 < 99 so this branch is unused
            end
        end
        # With the ceiling at 1.0 the same signals rail -- that is the knife-edge.
        @test solveAh(net, 1.5, 5000, "1.0").node_activities["A"] * 100 > 99.0
    end

    @testset "ceiling < 1 leaves the band, strong drive and knockout as they were" begin
        net = posloop()
        for hi in ("0.9", "0.95")
            @test isapprox(solveAh(net, 1.01, 2000, hi).node_activities["A"] * 100, 1.0106; rtol=1e-3)
            @test solveAh(net, 80.0, 500, hi).node_activities["A"] * 100 > 99.0
            @test solveAh(net, 0.0, 500, hi).node_activities["A"] * 100 < 1e-6
        end
        # Lower ceiling => less amplification above the band, monotone.
        r90 = solveAh(net, 1.5, 5000, "0.9").node_activities["A"]
        r95 = solveAh(net, 1.5, 5000, "0.95").node_activities["A"]
        @test r90 < r95
    end

    # ------------------------------------------------------------------
    # DS_LOOP_GAIN (specs/028): the pull is per LOOP, not per edge. Each
    # in-loop edge is read at fold^(g^(1/n)), n = its component's size, so a
    # cycle's total gain is g however long it is: A = U^(1/(1-g)).
    # ------------------------------------------------------------------
    solveG(net, U, it, g) = withenv("DS_LOOP_GAIN" => g) do
        DS.solve_steady_state(net, Dict("U" => (U, 1.0)), P(it))
    end
    function longloop()   # U -> A ; A -> r1 -> B -> r2 -> C -> r3 -> D -> r4 -> A (8 nodes)
        ids = ("U", "A", "r1", "B", "r2", "C", "r3", "D", "r4")
        nodes = Dict(id => N(id) for id in ids)
        edges = [E("U", "A", true, true, "input"),
                 E("A", "r1", true, true, "input"), E("r1", "B", true, true, "output"),
                 E("B", "r2", true, true, "input"), E("r2", "C", true, true, "output"),
                 E("C", "r3", true, true, "input"), E("r3", "D", true, true, "output"),
                 E("D", "r4", true, true, "input"), E("r4", "A", true, true, "output")]
        DS.ReactionNetwork(nodes, edges, Dict{String,DS.SetExpansionMapping}())
    end

    @testset "DS_LOOP_GAIN guard rails" begin
        for bad in ("0", "-0.1", "1.5", "abc", "nan")
            withenv("DS_LOOP_GAIN" => bad) do
                @test_throws ArgumentError DS.resolve_reaction_eval_config()
            end
        end
        withenv("DS_LOOP_GAIN" => "0.99", "DS_LOOP_ELASTICITY" => "0.9") do
            @test_throws ArgumentError DS.resolve_reaction_eval_config()   # one pull, not two
        end
    end

    @testset "DS_LOOP_GAIN = 1 is byte-identical to unset" begin
        for U in (1.01, 0.9, 80.0)
            ref = solveG(posloop(), U, 500, nothing)
            got = solveG(posloop(), U, 500, "1.0")
            @test got.node_activities == ref.node_activities
        end
    end

    @testset "a minimal loop gain turns the knife-edge into finite amplification" begin
        # Each trip round the cycle closes a fraction (1 - g) of the gap, so a
        # near-1 gain converges slowly: at g = 0.99 A reads 1.76 after 500
        # sweeps and 2.67 after 5,000. The root is unique and budget-independent
        # once converged -- but a minimal gain needs a large budget, which a
        # catalog arm must check (its convergence is reported).
        for (U, g, it) in ((1.01, "0.99", 5000), (1.1, "0.9", 500), (0.95, "0.95", 1000))
            want = U^(1 / (1 - parse(Float64, g)))
            a = solveG(posloop(), U, it, g).node_activities["A"] * 100
            a4 = solveG(posloop(), U, 4 * it, g).node_activities["A"] * 100
            @test isapprox(a4, want; rtol = 0.02)          # a unique root ...
            @test isapprox(a, a4; rtol = 0.02)             # ... the budget no longer moves
        end
        # Up stays up and grows; down stays down and deepens (Adam's intent).
        @test solveG(posloop(), 1.05, 5000, "0.95").node_activities["A"] * 100 > 1.05
        @test solveG(posloop(), 0.95, 5000, "0.95").node_activities["A"] * 100 < 0.95
    end

    @testset "the pull is per loop: a longer loop amplifies by the same factor" begin
        a4 = solveG(posloop(), 1.05, 5000, "0.9").node_activities["A"] * 100
        a8 = solveG(longloop(), 1.05, 5000, "0.9").node_activities["A"] * 100
        @test isapprox(a4, a8; rtol = 0.02)
        @test isapprox(a4, 1.05^10; rtol = 0.02)
    end

    @testset "DS_LOOP_GAIN leaves acyclic networks untouched" begin
        for U in (0.0, 0.5, 2.0, 50.0)
            obs = Dict("S" => (U, 1.0))
            a = withenv("DS_LOOP_GAIN" => nothing) do
                DS.solve_steady_state(chain(), obs, P(200))
            end
            b = withenv("DS_LOOP_GAIN" => "0.9") do
                DS.solve_steady_state(chain(), obs, P(200))
            end
            @test a.node_activities == b.node_activities
        end
    end

end

