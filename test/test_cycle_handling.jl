# Cycle handling — specs/008-cycle-handling
#
# These fixtures exist because the loop work was twice built on measurements
# that did not reproduce. They pin the behaviour that IS verified, and mark the
# one defect that is real with @test_broken so a fix announces itself instead of
# being silently absorbed.
#
# Baseline is 0.01 internal / 1.0 UI throughout.

using Test
using DeltaSignal
const DS = DeltaSignal

const LOOP_NODES = ("A", "B", "r1", "r2")

"""
Two-reaction cycle with two external supplies:
    U1 -> A ;  U2 -> B ;  A -> r1 -> B ;  B -> r2 -> A ;  B -> R
`and_join` selects AND vs OR aggregation at the points where each supply meets
the cycle's own recycled product — the distinction the whole feature turns on.
"""
function cycle_fixture(; and_join::Bool)
    N(id) = DS.NetworkNode(id, "R-HSA-" * id, "unknown", nothing, id, 0.01)
    nodes = Dict(id => N(id) for id in ("U1", "U2", "A", "B", "r1", "r2", "R"))
    E(s, t, and, pos, ty) = DS.LogicNetworkEdge(s, t, and, pos, 1.0, ty)
    j = and_join
    edges = [
        E("U1", "A", j, true, "input"),   E("U2", "B", j, true, "input"),
        E("A", "r1", true, true, "input"), E("r1", "B", j, true, "output"),
        E("B", "r2", true, true, "input"), E("r2", "A", j, true, "output"),
        E("B", "R", false, true, "input"),
    ]
    DS.ReactionNetwork(nodes, edges, Dict{String,DS.SetExpansionMapping}())
end

solve_ui(net, obs) = begin
    r = DS.solve_steady_state(net, obs,
            DS.SteadyStateParams(1.0, 0.1, 500, 1e-8, "penalty"))
    (Dict(k => v * 100.0 for (k, v) in r.node_activities), r)
end

ko(u) = Dict(u => (0.0, 1.0))

@testset "cycle handling" begin

    @testset "OR-joined cycle satisfies US1" begin
        net = cycle_fixture(and_join=false)

        v, r = solve_ui(net, Dict{String,Tuple{Float64,Float64}}())
        @test r.converged
        @test isapprox(v["A"], 1.0; atol=1e-6)   # unperturbed sits at baseline
        @test isapprox(v["B"], 1.0; atol=1e-6)

        # Scenario 1: one supply out, the survivor holds the cycle up.
        v, r = solve_ui(net, ko("U1"))
        @test r.converged
        @test v["A"] > 0.1    # NOT collapsed: two orders above the dead state
        @test v["B"] > 0.1
        @test v["A"] < 1.0    # but genuinely reduced

        # Scenario 2: raising the survivor scales the cycle with it.
        _, _ = solve_ui(net, ko("U1"))
        v10, _ = solve_ui(net, Dict("U1" => (0.0, 1.0), "U2" => (10.0, 1.0)))
        v50, _ = solve_ui(net, Dict("U1" => (0.0, 1.0), "U2" => (50.0, 1.0)))
        @test v50["B"] > v10["B"] > v["B"]

        # Scenario 3: inputs genuinely zero -> zero is still reachable.
        vz, _ = solve_ui(net, Dict("U1" => (0.0, 1.0), "U2" => (0.0, 1.0)))
        @test vz["A"] < 0.01
        @test vz["B"] < 0.01
    end

    @testset "AND-joined cycle collapses (the real defect)" begin
        # One external input 50x above baseline, yet the cycle sits far BELOW
        # baseline: the cycle's own recycled product is a required co-input, so
        # a knockout elsewhere cannot be compensated. Cycle-internal positive
        # edges are AND in generated networks, so this is the common case.
        net = cycle_fixture(and_join=true)
        v, r = solve_ui(net, Dict("U1" => (0.0, 1.0), "U2" => (50.0, 1.0)))
        @test r.converged                      # it converges — this is not a
                                               # convergence bug
        @test v["A"] < 0.01                    # current behaviour: collapsed
        @test_broken v["A"] > 1.0              # FR-001: should be held up
        @test_broken v["B"] > 1.0
    end

    @testset "the fixed point is unique — there is no wrong basin" begin
        # specs/008 was built on "all-zero is a spurious root the iteration
        # falls into". Iterating the identical map from seven starting states
        # reaches one root, so no anti-collapse guard can be justified as
        # basin selection.
        net = cycle_fixture(and_join=true)
        rxns = DS.convert_to_reaction_network(net)
        allnodes = collect(keys(net.nodes))
        u2i = Dict(u => i for (i, u) in enumerate(allnodes))
        bl = Dict(u => net.nodes[u].baseline for u in allnodes)
        idx, _, _ = DS.index_reactions(rxns, u2i, bl, Set{String}())
        cfg = DS.resolve_reaction_eval_config()

        function settle(start)
            x = fill(0.01, length(allnodes))
            x[u2i["U1"]] = 0.0
            x[u2i["U2"]] = 0.50
            for u in LOOP_NODES
                x[u2i[u]] = start
            end
            pinned = Set([u2i["U1"], u2i["U2"]])
            for _ in 1:20_000
                mx = 0.0
                for rxn in idx
                    rxn.target_idx in pinned && continue
                    f = DS.compute_reaction_output_vec(x, rxn; config=cfg)
                    d = abs(f - x[rxn.target_idx])
                    d > mx && (mx = d)
                    x[rxn.target_idx] = 0.5 * x[rxn.target_idx] + 0.5 * f
                end
                mx < 1e-14 && break
            end
            [x[u2i[u]] for u in LOOP_NODES]
        end

        reference = settle(0.0)
        for start in (0.01, 0.1, 0.3, 0.5, 0.8, 1.0)
            @test settle(start) ≈ reference atol = 1e-12
        end
    end

    @testset "convergence verdict agrees with the stopping rule" begin
        # FR-003 / FR-004. The 4.045534594765421e-6 recorded in specs/008 does
        # not reproduce: the residual tracks the tolerance and every
        # non-converged case clears with more iterations.
        net = cycle_fixture(and_join=false)
        obs = Dict("U2" => (80.0, 1.0))
        prev = Inf
        for tol in (1e-6, 1e-8, 1e-12)
            r = DS.solve_steady_state(net, obs,
                    DS.SteadyStateParams(1.0, 0.1, 2000, tol, "penalty"))
            @test r.converged
            @test r.final_residual < tol
            @test r.final_residual < prev    # tightening tol tightens the result
            prev = r.final_residual
        end
    end
    @testset "substrate edges are visible to cycle detection" begin
        # compute_reaction_output_vec reads x[substrate_indices] into the
        # availability factor L, so a substrate is a real dependency of the
        # target. It must therefore appear in the graph Tarjan runs on: an edge
        # missing there is a cycle the SCC solver cannot see, and a singleton
        # component is evaluated exactly once with no iteration to correct a
        # stale input.
        #
        # Nothing populates `substrate_uuids` today, so this is a latent
        # landmine rather than a live bug — which is exactly why it needs a
        # test rather than a measurement.
        params(ns) = DS.create_default_reaction_params(0, 0, ns)
        rxn(target, substrate) = DS.Reaction(
            target, String[], String[], String[], [substrate], String[],
            params(1), true, Bool[], Bool[])

        # Two reactions closing a cycle ONLY through substrate edges.
        reactions = [rxn("B", "A"), rxn("A", "B")]
        u2i = Dict("A" => 1, "B" => 2)
        bl = Dict("A" => 0.01, "B" => 0.01)
        _, comp_id, _ = DS.index_reactions(reactions, u2i, bl, Set{String}())

        @test comp_id[u2i["A"]] == comp_id[u2i["B"]]   # one component, not two
        @test count(==(comp_id[u2i["A"]]), comp_id) == 2
    end
end
