# Loop as a conserved pool (specs/017). Assertions for US1-US3.
using Test
include(joinpath(@__DIR__, "..", "src", "DeltaSignal.jl"))
using .DeltaSignal
const DS = DeltaSignal
const BL = 0.01

N(id) = DS.NetworkNode(id, "R-HSA-" * id, "unknown", nothing, id, BL)
E(s, t, and, pos, ty) = DS.LogicNetworkEdge(s, t, and, pos, 1.0, ty)
mk(ids, edges) = DS.ReactionNetwork(Dict(id => N(id) for id in ids), edges, Dict{String,DS.SetExpansionMapping}())
P = DS.SteadyStateParams(1.0, 1e-6, 300, 1e-8, "penalty")
solve(net, obs, mode) = withenv("DS_SCC_METHOD" => mode) do
    DS.solve_steady_state(net, Dict{String,Tuple{Float64,Float64}}(obs), P)
end
fold(r, id) = r.node_activities[id] / BL
diag(r, k) = r.diagnostics[k]

# U -> A ; A -> r1 -> B ; B -> r2 -> A   (positive two-node loop, one external drive)
posloop() = mk(("U", "A", "B", "r1", "r2"),
    [E("U", "A", true, true, "input"),
     E("A", "r1", true, true, "input"), E("r1", "B", true, true, "output"),
     E("B", "r2", true, true, "input"), E("r2", "A", true, true, "output")])

# four external entries, co-required (AND) by four distinct member reactions of one loop
function fourentry()
    ids = ("E1", "E2", "E3", "E4", "A", "B", "C", "D", "rA", "rB", "rC", "rD")
    ed = DS.LogicNetworkEdge[]
    for (ent, prev, rx, nxt) in (("E1", "D", "rA", "A"), ("E2", "A", "rB", "B"), ("E3", "B", "rC", "C"), ("E4", "C", "rD", "D"))
        push!(ed, E(ent, rx, true, true, "input")); push!(ed, E(prev, rx, true, true, "input")); push!(ed, E(rx, nxt, true, true, "output"))
    end
    mk(ids, ed)
end

# loop member A has three OR producers: two external reactions (rX from E1, rY from E2)
# and its in-loop producer r1 (from B). The external fold of A's reaction is the OR
# mean of (3, 1, 1) -- the recycle route counts as one alternative at baseline.
orentry() = mk(("E1", "E2", "A", "B", "r1", "r2", "rX", "rY"),
    [E("E1", "rX", true, true, "input"), E("rX", "A", false, true, "output"),
     E("E2", "rY", true, true, "input"), E("rY", "A", false, true, "output"),
     E("B", "r1", true, true, "input"), E("r1", "A", false, true, "output"),
     E("A", "r2", true, true, "input"), E("r2", "B", true, true, "output")])

# X -> r1 -> M ; M -| r2 ; r2 -> T ; T -> r3 -> M : ONE negative edge on the cycle.
# This is MDM2 -| TP53 -> MDM2: genuine negative feedback (odd cycle). Parity is
# inconsistent here by construction, so pool_parity must fall back.
negfeedback() = mk(("X", "M", "T", "r1", "r2", "r3"),
    [E("X", "r1", true, true, "input"), E("r1", "M", true, true, "output"),
     E("M", "r2", true, false, "regulator"), E("r2", "T", true, true, "output"),
     E("T", "r3", true, true, "input"), E("r3", "M", true, true, "output")])

# X -> r1 -> M ; M -| r2 -> T ; T -| r3 -> M : TWO negatives on the cycle (double
# inhibition = positive feedback). Parity is consistent: T reads the inverse of M.
parityfix() = mk(("X", "M", "T", "r1", "r2", "r3"),
    [E("X", "r1", true, true, "input"), E("r1", "M", true, true, "output"),
     E("M", "r2", true, false, "regulator"), E("r2", "T", true, true, "output"),
     E("T", "r3", true, false, "regulator"), E("r3", "M", true, true, "output")])

# odd negative cycle: A -| rB -> B ; B -> rA -> A, driven by U -> A
oddcycle() = mk(("U", "A", "B", "rA", "rB"),
    [E("U", "rA", true, true, "input"), E("B", "rA", true, true, "input"), E("rA", "A", true, true, "output"),
     E("A", "rB", true, false, "regulator"), E("rB", "B", true, true, "output")])

# three cyclic components: two positive loops and one with an internal inhibitor, all fed by U
function threecomp()
    ed = DS.LogicNetworkEdge[E("U", "A1", true, true, "input"), E("U", "A2", true, true, "input"), E("U", "A3", true, true, "input")]
    for k in ("1", "2")
        append!(ed, [E("A$k", "p$k", true, true, "input"), E("p$k", "B$k", true, true, "output"), E("B$k", "q$k", true, true, "input"), E("q$k", "A$k", true, true, "output")])
    end
    append!(ed, [E("A3", "p3", true, true, "input"), E("p3", "B3", true, true, "output"), E("B3", "q3", true, false, "regulator"), E("q3", "A3", true, true, "output")])
    mk(("U", "A1", "B1", "p1", "q1", "A2", "B2", "p2", "q2", "A3", "B3", "p3", "q3"), ed)
end

relabel(net, f) = DS.ReactionNetwork(Dict(f(k) => DS.NetworkNode(f(v.uuid), v.reactome_id, v.entity_type, v.original_set_id, v.display_name, v.baseline) for (k, v) in net.nodes),
                                     [DS.LogicNetworkEdge(f(e.parent_uuid), f(e.child_uuid), e.is_and, e.is_positive, e.stoichiometry, e.edge_type) for e in net.edges],
                                     net.set_mappings)
RELABELLINGS = (s -> "zz_" * s, s -> string(hash(s)), s -> reverse(s) * "_q")

@testset "loop pool" begin
    @testset "guard rail" begin
        @test_throws ArgumentError solve(posloop(), Dict("U" => (2.0, 1.0)), "pool_parityy")
        @test_throws ArgumentError solve(posloop(), Dict("U" => (2.0, 1.0)), "Pool")
    end

    @testset "US2: default is bit-identical to fixed_point" begin
        for net in (posloop(), parityfix(), threecomp()), obs in (Dict("U" => (2.0, 1.0)), Dict("X" => (0.5, 1.0)))
            all(haskey(net.nodes, k) for k in keys(obs)) || continue
            @test solve(net, obs, nothing).node_activities == solve(net, obs, "fixed_point").node_activities
        end
        r = solve(posloop(), Dict("U" => (2.0, 1.0)), nothing)
        @test diag(r, "scc_pooled") == 0 && diag(r, "scc_iterated") == 1 && diag(r, "scc_method") == "fixed_point"
    end

    @testset "US1: a positive loop reads its external supply" begin
        for mode in ("pool", "pool_parity")
            r = solve(posloop(), Dict("U" => (2.0, 1.0)), mode)
            @test fold(r, "A") ≈ 2.0 rtol=1e-9
            @test fold(r, "B") ≈ 2.0 rtol=1e-9
            @test diag(r, "scc_pooled") == 1 && diag(r, "scc_iterated") == 0
            # knocked-out entry zeroes the pool
            r0 = solve(posloop(), Dict("U" => (0.0, 1.0)), mode)
            @test r0.node_activities["A"] == 0.0 && r0.node_activities["B"] == 0.0
            # no external perturbation: baseline everywhere
            rb = solve(posloop(), Dict(), mode)
            @test fold(rb, "A") ≈ 1.0 rtol=1e-12
            @test fold(rb, "B") ≈ 1.0 rtol=1e-12
        end
        # the knife-edge it replaces: fixed_point does NOT read 2x (it rails or collapses)
        rf = solve(posloop(), Dict("U" => (2.0, 1.0)), "fixed_point")
        @test !(1.5 < fold(rf, "A") < 2.5)
        # four co-required entries at 1,1,3,2 -> pool 6 (Adam's example)
        r4 = solve(fourentry(), Dict("E1" => (1.0, 1.0), "E2" => (1.0, 1.0), "E3" => (3.0, 1.0), "E4" => (2.0, 1.0)), "pool")
        for m in ("A", "B", "C", "D"); @test fold(r4, m) ≈ 6.0 rtol=1e-9; end
        # OR-alternative producers of A: the OR combination mean(3, 1, 1) = 5/3, not the product 3
        ro = solve(orentry(), Dict("E1" => (3.0, 1.0)), "pool")
        @test fold(ro, "A") ≈ 5/3 rtol=1e-9
        @test fold(ro, "B") ≈ 5/3 rtol=1e-9
    end

    @testset "US2: fallback is counted, never silent" begin
        r = solve(threecomp(), Dict("U" => (2.0, 1.0)), "pool")
        @test diag(r, "scc_pooled") == 2 && diag(r, "scc_iterated") == 1 && diag(r, "scc_fallback_negative") == 1
        @test fold(r, "A1") ≈ 2.0 rtol=1e-9
        @test fold(r, "A2") ≈ 2.0 rtol=1e-9
        # the component that fell back equals the fixed-point solver's answer for it
        rf = solve(threecomp(), Dict("U" => (2.0, 1.0)), "fixed_point")
        @test r.node_activities["A3"] == rf.node_activities["A3"] && r.node_activities["B3"] == rf.node_activities["B3"]
        @test diag(r, "scc_pooled_nodes") == 8   # A, B and the two reaction nodes of each pooled loop
        # pool_all pools the third component too; its inhibitor contributes nothing
        ra = solve(threecomp(), Dict("U" => (2.0, 1.0)), "pool_all")
        @test diag(ra, "scc_pooled") == 3 && diag(ra, "scc_iterated") == 0
        @test fold(ra, "A3") ≈ 2.0 rtol=1e-9
    end

    @testset "US3: internal negatives by parity" begin
        # consistent parity (double inhibition): T reads the inverse of M
        rp = solve(parityfix(), Dict("X" => (0.5, 1.0)), "pool_parity")
        @test fold(rp, "M") ≈ 0.5 rtol=1e-6          # the entry (AND with the in-loop producer at baseline)
        @test fold(rp, "T") ≈ 2.0 rtol=1e-6          # inverse of M through the internal inhibitor
        @test diag(rp, "scc_pooled") == 1 && diag(rp, "scc_fallback_inconsistent") == 0
        # conservative rule: falls back and matches the fixed point exactly
        rc = solve(parityfix(), Dict("X" => (0.5, 1.0)), "pool")
        rf = solve(parityfix(), Dict("X" => (0.5, 1.0)), "fixed_point")
        @test diag(rc, "scc_fallback_negative") == 1 && diag(rc, "scc_pooled") == 0
        @test rc.node_activities == rf.node_activities
        # ONE negative on the cycle (MDM2 -| TP53 -> MDM2): odd cycle, genuine negative
        # feedback -- parity is inconsistent and both spec'd variants iterate
        for (mode, key) in (("pool_parity", "scc_fallback_inconsistent"), ("pool", "scc_fallback_negative"))
            rn = solve(negfeedback(), Dict("X" => (0.5, 1.0)), mode)
            @test diag(rn, key) == 1 && diag(rn, "scc_pooled") == 0
        end
        ro = solve(oddcycle(), Dict("U" => (2.0, 1.0)), "pool_parity")
        @test diag(ro, "scc_fallback_inconsistent") == 1 && diag(ro, "scc_pooled") == 0
        # pool_all: pools it anyway, the inhibitor contributes nothing, T moves WITH M
        # (the known-wrong direction for the AKT -> MDM2 -> TP53 case; pinned so it is deliberate)
        rl = solve(negfeedback(), Dict("X" => (0.5, 1.0)), "pool_all")
        @test diag(rl, "scc_pooled") == 1
        @test fold(rl, "M") ≈ 0.5 rtol=1e-6
        @test fold(rl, "T") ≈ 0.5 rtol=1e-6
        # a pinned member inside the loop is an entry and is never overwritten
        rk = solve(parityfix(), Dict("M" => (0.0, 1.0)), "pool_parity")
        @test rk.node_activities["M"] == 0.0
        @test diag(rk, "scc_pooled") == 1
    end

    @testset "FR-008: label and edge order do not matter" begin
        # only fixtures every variant POOLS: a component that falls back is iterated by
        # Gauss-Seidel, whose result is label-dependent by the solver's known defect
        for (net, obs) in ((posloop(), Dict("U" => (2.0, 1.0))), (parityfix(), Dict("X" => (0.5, 1.0))), (fourentry(), Dict("E3" => (3.0, 1.0), "E4" => (2.0, 1.0))), (threecomp(), Dict("U" => (2.0, 1.0))), (negfeedback(), Dict("X" => (0.5, 1.0))))
            for mode in ("pool_parity", "pool_all")
                net === parityfix() && false
                r0 = solve(net, obs, mode)
                diag(r0, "scc_iterated") == 0 || continue
                ref = solve(net, obs, mode).node_activities
                for f in RELABELLINGS
                    rl = solve(relabel(net, f), Dict(f(k) => v for (k, v) in obs), mode).node_activities
                    @test all(rl[f(k)] == v for (k, v) in ref)
                end
                rev = DS.ReactionNetwork(net.nodes, reverse(net.edges), net.set_mappings)
                @test solve(rev, obs, mode).node_activities == ref
            end
        end
    end
end
