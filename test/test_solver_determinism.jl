# Solver determinism under node relabelling
#
# The solver's answer must be a function of the GRAPH, not of the names its
# nodes happen to carry. It was not: node indices come from
# `collect(keys(network.nodes))`, i.e. Julia `Dict` order, i.e. UUID hashing.
# That ordering reaches the Gauss-Seidel sweep inside a cyclic component, where
# each reaction reads values its neighbours have already updated this pass — so
# the update order selects WHICH fixed point a multi-root component lands in.
#
# Measured on the real catalog: renaming every UUID in the 92-pathway set (the
# relabelling verified isomorphic — identical stable-id canonical form for all
# 92) moved 14 of 23,908 curator predictions. Every one was in a pathway with a
# cycle; no acyclic pathway moved a single case.
#
# `DS_SCC_SWEEP=jacobi` evaluates a sweep against the state at its start and
# commits together, which removes the order dependence. These tests pin that:
# Jacobi is invariant to relabelling, and acyclic networks are invariant either
# way (their solve is exact and order cannot matter).
#
# SCOPE OF THESE FIXTURES, stated plainly. The cyclic fixture below has ONE
# fixed point, and on it Gauss-Seidel is not exactly invariant either -- it
# deviates across relabellings by 2e-9 to 4e-9, i.e. at the 1e-8 convergence
# tolerance, not at the level that could flip a prediction. So these tests
# demonstrate that Jacobi is EXACTLY order-free where Gauss-Seidel is only
# approximately so; they do NOT reproduce the catalog's prediction flips. Those
# need a component that does not converge inside the iteration budget, which a
# seven-node fixture will not do. Do not read a green run here as evidence that
# the catalog-scale behaviour is fixed.
#
# This does NOT fix the underlying multiple-fixed-point problem (the all-zero
# root, specs/004). It removes node naming as the thing that chooses between
# the roots.

using Test
using DeltaSignal
const DS = DeltaSignal

"""Relabel every node id via `f`, preserving the graph exactly."""
function relabel(net::DS.ReactionNetwork, f)
    nodes = Dict(f(k) => DS.NetworkNode(f(v.uuid), v.reactome_id, v.entity_type,
                                        v.original_set_id, v.display_name, v.baseline)
                 for (k, v) in net.nodes)
    # a colliding relabelling would silently merge nodes and compare a different graph
    length(nodes) == length(net.nodes) || error("relabelling collided: $(length(net.nodes)) -> $(length(nodes)) nodes")
    edges = [DS.LogicNetworkEdge(f(e.parent_uuid), f(e.child_uuid), e.is_and,
                                 e.is_positive, e.stoichiometry, e.edge_type)
             for e in net.edges]
    DS.ReactionNetwork(nodes, edges, net.set_mappings, net.cofactor_stids)
end

"""Two-reaction cycle fed by two external supplies, plus a readout off B."""
function cyclic_net()
    N(id) = DS.NetworkNode(id, "R-HSA-" * id, "unknown", nothing, id, 0.01)
    nodes = Dict(id => N(id) for id in ("U1","U2","A","B","r1","r2","R"))
    E(s,t,and,pos,ty) = DS.LogicNetworkEdge(s,t,and,pos,1.0,ty)
    edges = [
        E("U1","A",false,true,"input"),  E("U2","B",false,true,"input"),
        E("A","r1",true,true,"input"),   E("r1","B",false,true,"output"),
        E("B","r2",true,true,"input"),   E("r2","A",false,true,"output"),
        E("B","R",false,true,"input"),
    ]
    DS.ReactionNetwork(nodes, edges, Dict{String,DS.SetExpansionMapping}())
end

"""Straight chain, no cycle anywhere."""
function acyclic_net()
    N(id) = DS.NetworkNode(id, "R-HSA-" * id, "unknown", nothing, id, 0.01)
    nodes = Dict(id => N(id) for id in ("S","r1","M","r2","T"))
    E(s,t,and,pos,ty) = DS.LogicNetworkEdge(s,t,and,pos,1.0,ty)
    edges = [E("S","r1",true,true,"input"), E("r1","M",false,true,"output"),
             E("M","r2",true,true,"input"), E("r2","T",false,true,"output")]
    DS.ReactionNetwork(nodes, edges, Dict{String,DS.SetExpansionMapping}())
end

params() = DS.SteadyStateParams(1.0, 0.1, 500, 1e-8, "penalty")

"""Solve and return the activities keyed by the ORIGINAL (pre-relabel) name."""
function solve_by_original(net, obs, f)
    r = DS.solve_steady_state(relabel(net, f),
                              Dict(f(k) => v for (k, v) in obs), params())
    Dict(k => r.node_activities[f(k)] for k in keys(net.nodes))
end

# A handful of structurally different relabellings. `pad` and `rev` change the
# hash of every id; `swap` additionally permutes which name lands on which node
# among the cycle's own members, which is what reorders the sweep.
const RELABELLINGS = Dict(
    "identity" => identity,
    "prefix"   => (s -> "zz_" * s),
    "suffix"   => (s -> s * "_9f3c"),
    "reversed" => (s -> reverse(s)),
    "uuidish"  => (s -> string("f47ac10b-58cc-4372-a567-0e02b2c3d4", lpad(sum(codeunits(s)) % 100, 2, '0'))),
)

withenv_sweep(f, mode) = withenv(f, "DS_SCC_SWEEP" => mode)

@testset "solver determinism" begin

    @testset "DS_SCC_SWEEP validation" begin
        net, obs = acyclic_net(), Dict("S" => (50.0, 1.0))
        @test_throws ArgumentError withenv_sweep(
            () -> DS.solve_steady_state(net, obs, params()), "gauss-seidel")
        @test_throws ArgumentError withenv_sweep(
            () -> DS.solve_steady_state(net, obs, params()), "")
        for good in ("gauss_seidel", "jacobi")
            @test withenv_sweep(
                () -> DS.solve_steady_state(net, obs, params()), good) !== nothing
        end
    end

    @testset "acyclic networks are relabel-invariant in both sweeps" begin
        net = acyclic_net()
        for mode in ("gauss_seidel", "jacobi"), obsval in (0.0, 1.0, 50.0, 100.0)
            obs = Dict("S" => (obsval, 1.0))
            ref = withenv_sweep(() -> solve_by_original(net, obs, identity), mode)
            for (nm, f) in RELABELLINGS
                got = withenv_sweep(() -> solve_by_original(net, obs, f), mode)
                @test all(isapprox(got[k], ref[k]; atol=0, rtol=0) for k in keys(ref))
            end
        end
    end

    @testset "jacobi is relabel-invariant on a cyclic network" begin
        net = cyclic_net()
        for obsval in (0.0, 1.0, 50.0, 100.0)
            obs = Dict("U1" => (obsval, 1.0))
            ref = withenv_sweep(() -> solve_by_original(net, obs, identity), "jacobi")
            for (nm, f) in RELABELLINGS
                got = withenv_sweep(() -> solve_by_original(net, obs, f), "jacobi")
                @test all(isapprox(got[k], ref[k]; atol=0, rtol=0) for k in keys(ref))
            end
        end
    end

    @testset "jacobi reaches the same fixed point as gauss_seidel here" begin
        # This fixture's component has a single root, so the two schemes must
        # agree on it. The point is to guard against Jacobi being order-free by
        # being WRONG -- converging somewhere else, or not at all -- rather than
        # by being correct. Agreement to 1e-6 is the check; exact agreement is
        # not expected, since the two take different paths to the same root.
        net = cyclic_net()
        for obsval in (0.0, 1.0, 50.0, 100.0)
            obs = Dict("U1" => (obsval, 1.0))
            gs = withenv_sweep(() -> solve_by_original(net, obs, identity), "gauss_seidel")
            ja = withenv_sweep(() -> solve_by_original(net, obs, identity), "jacobi")
            @test all(isapprox(ja[k], gs[k]; atol=1e-6) for k in keys(gs))
        end
    end

    @testset "observations stay pinned under both sweeps" begin
        net = cyclic_net()
        # Observations arrive on the UI 0-100 scale; activities come back
        # internal, so the pinned value is v/100.
        for mode in ("gauss_seidel", "jacobi"), v in (0.0, 25.0, 100.0)
            r = withenv_sweep(() -> DS.solve_steady_state(
                    net, Dict("U1" => (v, 1.0), "U2" => (v, 1.0)), params()), mode)
            @test isapprox(r.node_activities["U1"], v / 100.0; atol=1e-9)
            @test isapprox(r.node_activities["U2"], v / 100.0; atol=1e-9)
        end
    end
end
