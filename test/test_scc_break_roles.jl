# Derived-edge recycling closures do not carry loop signal (specs/018).
using Test
include(joinpath(@__DIR__, "..", "src", "DeltaSignal.jl"))
using .DeltaSignal
const DS = DeltaSignal
const BL = 0.01

N(id) = DS.NetworkNode(id, "R-HSA-" * id, "unknown", nothing, id, BL)
E(s, t, and, pos, ty) = DS.LogicNetworkEdge(s, t, and, pos, 1.0, ty)
mk(ids, edges) = DS.ReactionNetwork(Dict(id => N(id) for id in ids), edges, Dict{String,DS.SetExpansionMapping}())
P = DS.SteadyStateParams(1.0, 1e-6, 400, 1e-8, "penalty")
solve(net, obs, roles) = withenv("DS_SCC_BREAK_ROLES" => roles) do
    DS.solve_steady_state(net, Dict{String,Tuple{Float64,Float64}}(obs), P)
end
fold(r, id) = r.node_activities[id] / BL
diag(r, k) = r.diagnostics[k]

# Two reaction cycles welded by ONE assembly edge:
#   cycle 1: U -> A ; A -> r1 -> B ; B -> r2 -> A            (B is a member released by cycle 1)
#   cycle 2: C -> r3 -> D ; D -> r4 -> C                      (C is a complex)
#   weld:    B --(assembly)--> C   and   D --(input)--> r1    (so C's cycle reaches back into cycle 1)
# One SCC {A,B,r1,r2,C,D,r3,r4} today; with assembly closures broken, two.
function welded()
    ids = ("U", "A", "B", "r1", "r2", "C", "D", "r3", "r4", "M")
    mk(ids, [E("U", "A", true, true, "input"),
             E("A", "r1", true, true, "input"), E("r1", "B", true, true, "output"),
             E("B", "r2", true, true, "input"), E("r2", "A", true, true, "output"),
             E("C", "r3", true, true, "input"), E("r3", "D", true, true, "output"),
             E("D", "r4", true, true, "input"), E("r4", "C", true, true, "output"),
             E("B", "C", true, true, "assembly"), E("D", "r1", true, true, "input"),
             E("M", "C", true, true, "assembly")])            # M: an external member (feed-forward, not a closure)
end
# depletion closure: complex X drains its subunit S inside one cycle; Z drives S from outside
depweld() = mk(("Z", "S", "rS", "X", "rX", "Y", "rY"),
    [E("Z", "S", true, true, "input"),
     E("S", "rX", true, true, "input"), E("rX", "X", true, true, "output"),
     E("X", "rY", true, true, "input"), E("rY", "Y", true, true, "output"),
     E("Y", "rS", true, true, "input"), E("rS", "S", true, true, "output"),
     E("X", "S", true, false, "depletion")])
relabel(net, f) = DS.ReactionNetwork(Dict(f(k) => DS.NetworkNode(f(v.uuid), v.reactome_id, v.entity_type, v.original_set_id, v.display_name, v.baseline) for (k, v) in net.nodes),
                                     [DS.LogicNetworkEdge(f(e.parent_uuid), f(e.child_uuid), e.is_and, e.is_positive, e.stoichiometry, e.edge_type) for e in net.edges], net.set_mappings)

@testset "SCC break roles" begin
    @testset "guard rail" begin
        @test_throws ArgumentError solve(welded(), Dict("U" => (2.0, 1.0)), "assembly,depletio")
        @test_throws ArgumentError solve(welded(), Dict("U" => (2.0, 1.0)), "Assembly")
    end

    @testset "US3: unset is bit-identical; the legacy catalyst knob is untouched" begin
        for net in (welded(), depweld()), obs in (Dict("U" => (2.0, 1.0)), Dict("Z" => (0.5, 1.0)))
            all(haskey(net.nodes, k) for k in keys(obs)) || continue
            @test solve(net, obs, nothing).node_activities == solve(net, obs, "").node_activities
        end
        r = solve(welded(), Dict("U" => (2.0, 1.0)), nothing)
        @test diag(r, "scc_closures_assembly") == 0 && diag(r, "scc_cyclic_before") == diag(r, "scc_cyclic_after")
    end

    @testset "US1: a welded component falls apart, feed-forward limiting survives" begin
        r = solve(welded(), Dict("U" => (2.0, 1.0)), "assembly")
        @test diag(r, "scc_cyclic_before") == 1
        @test diag(r, "scc_cyclic_after") == 2
        @test diag(r, "scc_closures_assembly") == 1      # B -> C only; M -> C is not a closure
        @test diag(r, "scc_closures_catalyst") == 0 && diag(r, "scc_closures_depletion") == 0
        # an external member knocked out still lowers the complex (feed-forward assembly kept)
        rk = solve(welded(), Dict("M" => (0.0, 1.0)), "assembly")
        @test rk.node_activities["C"] < 1e-9
        # the rule leaves a component intact when its closures are of another role
        r2 = solve(welded(), Dict("U" => (2.0, 1.0)), "catalyst")
        @test diag(r2, "scc_cyclic_after") == 1 && diag(r2, "scc_closures_catalyst") == 0
    end

    @testset "US1: depletion closure reads the entry value, external drive still suppresses" begin
        r = solve(depweld(), Dict("Z" => (2.0, 1.0)), "depletion")
        @test diag(r, "scc_closures_depletion") == 1
        @test diag(r, "scc_cyclic_after") == 1        # the reaction cycle S -> X -> Y -> S remains
        # a depletion edge from OUTSIDE the component is untouched: raise X's cycle-external twin
        ext = mk(("Q", "S", "W", "rW"), [E("Q", "W", true, true, "input"), E("W", "rW", true, true, "input"), E("rW", "S", true, true, "output"), E("W", "S", true, false, "depletion")])
        a = solve(ext, Dict("Q" => (20.0, 1.0)), "depletion").node_activities["S"]
        b = solve(ext, Dict("Q" => (20.0, 1.0)), nothing).node_activities["S"]
        @test a == b
    end

    @testset "FR-006: closure marking is label- and edge-order-free; remnants iterate as before" begin
        # The closure SET and the recomputed component census must be identical under
        # relabelling and edge reversal (exact). Activities of components that still
        # iterate carry the solver's known Gauss-Seidel label drift (~1e-8, specs/013)
        # with or without roles, so they are compared to solver tolerance, not bitwise.
        for (net, obs) in ((welded(), Dict("U" => (2.0, 1.0))), (depweld(), Dict("Z" => (0.5, 1.0))))
            for roles in ("assembly", "assembly,depletion", "catalyst,assembly,depletion")
                r0 = solve(net, obs, roles)
                census(r) = (diag(r, "scc_closures_assembly"), diag(r, "scc_closures_depletion"), diag(r, "scc_closures_catalyst"), diag(r, "scc_cyclic_before"), diag(r, "scc_cyclic_after"), diag(r, "scc_largest_after"))
                for f in (s -> "zz_" * s, s -> string(hash(s)))
                    rl = solve(relabel(net, f), Dict(f(k) => v for (k, v) in obs), roles)
                    @test census(rl) == census(r0)
                    @test all(abs(rl.node_activities[f(k)] - v) <= 1e-7 for (k, v) in r0.node_activities)
                end
                rev = DS.ReactionNetwork(net.nodes, reverse(net.edges), net.set_mappings)
                rr = solve(rev, obs, roles)
                @test census(rr) == census(r0)
                @test all(abs(rr.node_activities[k] - v) <= 1e-7 for (k, v) in r0.node_activities)
            end
        end
    end

    @testset "review fixes: closures read the entry state on every evaluation path; self-loops are not closures" begin
        # a closure activator into a POOLED component must carry an upstream perturbation
        ko = Dict{String,Tuple{Float64,Float64}}("M" => (0.0, 1.0))
        r_fp = withenv("DS_SCC_BREAK_ROLES" => "assembly") do
            DS.solve_steady_state(welded(), ko, P)
        end
        r_pool = withenv("DS_SCC_BREAK_ROLES" => "assembly", "DS_SCC_METHOD" => "pool_all") do
            DS.solve_steady_state(welded(), ko, P)
        end
        @test r_fp.node_activities["C"] < 1e-9
        @test r_pool.node_activities["C"] < 1e-9
        # the flat solver honours closures the same way
        r_flat = withenv("DS_SCC_BREAK_ROLES" => "assembly", "DS_SCC_SOLVE" => "0") do
            DS.solve_steady_state(welded(), ko, P)
        end
        @test r_flat.node_activities["C"] < 1e-9
        # influence of the limiting member is not zeroed by a closure clamped at fold 1
        r = solve(welded(), Dict("U" => (2.0, 1.0)), "assembly")
        infl = withenv("DS_SCC_BREAK_ROLES" => "assembly") do
            DS.compute_influence_scores(r, DS.convert_to_reaction_network(welded()))
        end
        @test get(infl, "M", 0.0) > 0.0
        # a self-catalysing reaction is not a closure
        # a DIRECT self-edge (node catalyses its own production); the two-node A <-> rA form is a real cycle
        selfcat = mk(("U", "A"), [E("U", "A", true, true, "input"), E("A", "A", true, true, "catalyst")])
        rs = solve(selfcat, Dict("U" => (2.0, 1.0)), "catalyst")
        @test diag(rs, "scc_closures_catalyst") == 0
        @test diag(rs, "scc_cyclic_before") == diag(rs, "scc_cyclic_after")
    end
end
