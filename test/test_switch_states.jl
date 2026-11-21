# specs/154: switch detection — a second stable fixed point per cyclic
# component, served beside the canonical answer, never instead of it.
#
# The bistable fixture COMPUTES its two states (nothing pins the loop): driver
# P feeds an AND reaction rA together with loop partner B; rA makes A; rB reads
# A and makes B. Under hill_sat AND the fold product gives fold_A = fold_P *
# fold_B, fold_B = fold_A, so with P pinned at 2x the loop's only
# self-consistent states are the all-zero state and the cap — exactly the
# knife-edge class (specs/014) behind the recorded rails. A readout D inhibited
# by B tests the per-readout "high"/"low" naming (D is LOW when the loop is
# high). The monostable control is the same shape with the loop edge
# inhibitory (A = P / B has one positive fixed point), where the re-solves
# must come home and claim nothing.
using Test
using DeltaSignal
const DS = DeltaSignal

function with_env(f, pairs...)
    old = Dict(k => get(ENV, k, nothing) for (k, _) in pairs)
    for (k, v) in pairs
        v === nothing ? delete!(ENV, k) : (ENV[k] = v)
    end
    try
        f()
    finally
        for (k, v) in old
            v === nothing ? delete!(ENV, k) : (ENV[k] = v)
        end
    end
end

const BL = 0.01
e(s, t, et = "input") = DS.LogicNetworkEdge(s, t, et != "output", et != "depletion" && et != "inhibitor", 1.0,
                                            et == "inhibitor" ? "regulator_neg" : et)

# loop_sign = :pos — B feeds rA as a second AND input (bistable with P > 1x).
# loop_sign = :neg — B inhibits rA (one fixed point).
function loopnet(; loop_sign = :pos, with_readout = true)
    names = ["P", "A", "B", "rA", "rB"]
    with_readout && append!(names, ["C0", "D", "rD"])
    nodes = Dict(u => DS.NetworkNode(u, "R-$u", "protein", nothing, u, BL) for u in names)
    edges = [e("P", "rA"), e("rA", "A", "output"),
             e("A", "rB"), e("rB", "B", "output"),
             loop_sign == :pos ? e("B", "rA") : e("B", "rA", "inhibitor")]
    if with_readout
        append!(edges, [e("C0", "rD"), e("B", "rD", "inhibitor"), e("rD", "D", "output")])
    end
    DS.ReactionNetwork(nodes, edges, Dict{String, DS.SetExpansionMapping}(), Set{String}(),
                       Dict{String, Set{String}}(), nothing, nothing, nothing, nothing, nothing)
end

# Review 2026-11-20 (B1): the monostable loop plus a member R that no reaction
# writes. R is an AND input of rA and joins the cycle only through a structural
# `recycle` edge B -> R (specs/081), exactly as RBL2's starting copy does in
# Mitotic G1. The sweep never recomputes R, so a seed there would survive.
function recyclenet()
    names = ["P", "A", "B", "rA", "rB", "R"]
    nodes = Dict(u => DS.NetworkNode(u, "R-$u", "protein", nothing, u, BL) for u in names)
    edges = [e("P", "rA"), e("R", "rA"), e("rA", "A", "output"),
             e("A", "rB"), e("rB", "B", "output"), e("B", "rA", "inhibitor"),
             DS.LogicNetworkEdge("B", "R", true, true, 1.0, "recycle")]
    DS.ReactionNetwork(nodes, edges, Dict{String, DS.SetExpansionMapping}(), Set{String}(),
                       Dict{String, Set{String}}(), nothing, nothing, nothing, nothing, nothing)
end

# Review 2026-11-20 (B2): the bistable loop beside an UNRELATED autocatalytic
# loop (Q drives rN; its product N catalyses rN). Nothing connects the two, so
# a flip of the switch cannot move N. Under DS_SELF_FED_MODE=catalyst_loop
# (specs/109) the N -> rN catalyst edge is held at the component's ENTRY value,
# so canonical N reads Q's 4x. Re-sweeping that loop with the edge read live
# (what the old post-pass did to every lower-numbered component) rails it to
# the cap, and the old code served the move as part of the switch's alternate.
# Component numbers follow node-name hashing, so `prefix` renames the
# unrelated loop to put it on both sides of the switch in topological order
# (with "M" it lands below the switch; the old code failed there).
function twoloops(prefix)
    net = loopnet(with_readout = false)
    q, n_, rn = (prefix * x for x in ("Q", "N", "rN"))
    for u in (q, n_, rn)
        net.nodes[u] = DS.NetworkNode(u, "R-$u", "protein", nothing, u, BL)
    end
    append!(net.edges, [e(q, rn), e(rn, n_, "output"), e(n_, rn, "catalyst")])
    return net, (q, n_, rn)
end

solve(net, obs; iters = 500) =
    DS.solve_steady_state(net, obs, DS.SteadyStateParams(1.0, 0.1, iters, 1e-6, "penalty"))

sw_comps(r) = get(r.diagnostics, "switch_components", Vector{Dict{String, Any}}())

@testset "switch states (specs/154)" begin

@testset "env flags validate and default" begin
    with_env("DS_SWITCH_DETECT" => nothing) do
        @test DS._switch_detect_on()
    end
    with_env("DS_SWITCH_DETECT" => "off") do
        @test !DS._switch_detect_on()
    end
    with_env("DS_SWITCH_DETECT" => "ON") do
        @test_throws ArgumentError DS._switch_detect_on()
    end
    with_env("DS_SWITCH_MAX_SWEEPS" => nothing) do
        @test DS._switch_max_sweeps() == 250
    end
    with_env("DS_SWITCH_MAX_SWEEPS" => "0") do
        @test_throws ArgumentError DS._switch_max_sweeps()
    end
    with_env("DS_SWITCH_FLOOR" => nothing) do
        @test DS._switch_floor() == 1e-3
    end
    with_env("DS_SWITCH_FLOOR" => "-1") do
        @test_throws ArgumentError DS._switch_floor()
    end
end

@testset "the bistable loop: two states found, canonical untouched (US1)" begin
    r = solve(loopnet(), Dict("P" => (2.0, 1.0)))
    @test r.converged
    # canonical: the loop rails high (gain 2 from baseline grows to the cap)
    @test r.node_activities["A"] > 0.9
    @test r.node_activities["B"] > 0.9
    comps = sw_comps(r)
    @test r.diagnostics["switch_detect"] == "on"
    @test r.diagnostics["switch_components_found"] == 1
    @test length(comps) == 1
    alt = comps[1]["alt_activities"]
    # the alternate is the low state, for the loop members
    @test haskey(alt, "A") && alt["A"] < 1e-3
    @test haskey(alt, "B") && alt["B"] < 1e-3
    # and the flip PROPAGATES: D (inhibited by B) differs between the states
    @test haskey(alt, "D")
    @test alt["D"] > r.node_activities["D"]   # loop off => D de-repressed
    # the readout rule names states by the readout's OWN value
    fold_D = r.node_activities["D"] / BL
    sw, st, osf = DS.readout_switch(fold_D, Union{Nothing, Float64}[alt["D"]], BL)
    @test sw === true
    @test st == "low"            # canonical D is the lower value
    @test osf ≈ alt["D"] / BL
end

@testset "detection changes nothing canonical (FR-001, fixture level)" begin
    obs = Dict("P" => (2.0, 1.0))
    r_on  = solve(loopnet(), obs)
    r_off = with_env("DS_SWITCH_DETECT" => "off") do
        solve(loopnet(), obs)
    end
    @test r_off.diagnostics["switch_detect"] == "off"
    @test isempty(sw_comps(r_off))
    for (k, v) in r_off.node_activities
        @test r_on.node_activities[k] == v
    end
    @test r_on.converged == r_off.converged
    @test r_on.final_residual == r_off.final_residual
end

@testset "monostable control claims nothing (US2)" begin
    r = solve(loopnet(loop_sign = :neg), Dict("P" => (2.0, 1.0)))
    @test r.converged
    @test r.diagnostics["switch_components_found"] == 0
    @test isempty(sw_comps(r))
end

@testset "a pin inside the component removes the second state (edge case)" begin
    # B pinned low: A = 2 * 0 = 0 is the only consistent state; nothing to find.
    r = solve(loopnet(), Dict("P" => (2.0, 1.0), "B" => (0.0, 1.0)))
    @test r.converged
    @test r.diagnostics["switch_components_found"] == 0
end

@testset "the noise floor rejects near-identical points (FR-007)" begin
    # with the floor far above the whole scale, the distinct low state is
    # rejected as "the same point"
    r = with_env("DS_SWITCH_FLOOR" => "2.0") do
        solve(loopnet(), Dict("P" => (2.0, 1.0)))
    end
    @test r.converged
    @test r.diagnostics["switch_components_found"] == 0
end

@testset "honest absence under a starved sweep budget (FR-004/FR-008)" begin
    r = with_env("DS_SWITCH_MAX_SWEEPS" => "1") do
        solve(loopnet(), Dict("P" => (2.0, 1.0)))
    end
    @test r.converged
    # From the all-zero seed the bistable loop is already at its low fixed
    # point (A = P * 0 = 0) and settles in one sweep; from the all-one seed it
    # is at the cap. Neither basin is capped, and the low state is still found.
    @test r.diagnostics["switch_sweeps_capped"] == 0
    @test r.diagnostics["switch_components_found"] == 1
    # The monostable loop needs many sweeps from either seed: with one sweep
    # both basins are capped and counted, and nothing is served.
    r2 = with_env("DS_SWITCH_MAX_SWEEPS" => "1") do
        solve(loopnet(loop_sign = :neg), Dict("P" => (2.0, 1.0)))
    end
    @test r2.converged
    @test r2.diagnostics["switch_components_searched"] == 1
    @test r2.diagnostics["switch_sweeps_capped"] == 2
    @test r2.diagnostics["switch_components_found"] == 0
    @test isempty(sw_comps(r2))
end

@testset "a member no reaction writes is not seeded (review 2026-11-20, B1)" begin
    r = solve(recyclenet(), Dict("P" => (2.0, 1.0)))
    @test r.converged
    @test r.diagnostics["recycle_edges"] == 1
    # R has no producer: canonical keeps its baseline
    @test r.node_activities["R"] == BL
    # The loop is monostable; R seeded at 1.0 or 0.0 would have stood untouched
    # and the loop would have re-settled around it, a fabricated second state.
    @test r.diagnostics["switch_components_searched"] == 1
    @test r.diagnostics["switch_components_found"] == 0
    @test isempty(sw_comps(r))
end

@testset "an unrelated component is not part of the alternate (review 2026-11-20, B2)" begin
    for prefix in ("M", "a", "z", "Z", "0", "x_", "q")
        net, unrelated = twoloops(prefix)
        obs = Dict("P" => (2.0, 1.0), unrelated[1] => (4.0, 1.0))
        r = with_env("DS_SELF_FED_MODE" => "catalyst_loop") do
            solve(net, obs)
        end
        r_off = with_env("DS_SELF_FED_MODE" => "catalyst_loop", "DS_SWITCH_DETECT" => "off") do
            solve(net, obs)
        end
        @test r.converged
        @test r.diagnostics["self_fed_edges_held"] == 1
        # canonical: the held catalyst reads its entry (baseline), so N = Q = 4x
        @test isapprox(r.node_activities[unrelated[2]], 4 * BL; rtol = 1e-3)
        comps = sw_comps(r)
        @test r.diagnostics["switch_components_found"] == 1
        @test length(comps) == 1
        alt = comps[1]["alt_activities"]
        # the switch's own low state is served
        @test haskey(alt, "A") && alt["A"] < 1e-3
        @test haskey(alt, "B") && alt["B"] < 1e-3
        # nothing of the unrelated loop rides along
        for u in unrelated
            @test !haskey(alt, u)
        end
        # and the switch pass leaves the canonical answer as it was
        for (k, v) in r_off.node_activities
            @test r.node_activities[k] == v
        end
    end
end

@testset "a downstream component that does not settle caps the alternate (review 2026-11-20, B2)" begin
    # The switch's B drives a monostable negative loop X/Y downstream. With one
    # sweep, the switch's low basin settles at once (A = P * 0 = 0) but the
    # downstream loop, re-solved under B = 0, cannot: the alternate is capped
    # and counted, never served with the downstream loop's drifting values.
    function downnet()
        net = loopnet(with_readout = false)
        for u in ("X", "Y", "rX", "rY")
            net.nodes[u] = DS.NetworkNode(u, "R-$u", "protein", nothing, u, BL)
        end
        append!(net.edges, [e("B", "rX"), e("C1", "rX"), e("rX", "X", "output"), e("X", "rY"),
                            e("rY", "Y", "output"), e("Y", "rX", "inhibitor")])
        net.nodes["C1"] = DS.NetworkNode("C1", "R-C1", "protein", nothing, "C1", BL)
        return net
    end
    obs = Dict("P" => (2.0, 1.0), "C1" => (3.0, 1.0))
    r = solve(downnet(), obs)
    @test r.converged
    @test r.diagnostics["switch_components_found"] == 1     # full budget: found
    @test haskey(sw_comps(r)[1]["alt_activities"], "X")      # and propagated downstream
    r1 = with_env("DS_SWITCH_MAX_SWEEPS" => "1") do
        solve(downnet(), obs)
    end
    @test r1.converged
    @test r1.diagnostics["switch_components_found"] == 0
    @test r1.diagnostics["switch_sweeps_capped"] >= 1
    @test isempty(sw_comps(r1))
end

@testset "readout_switch: the call-cutoff rule (FR-003)" begin
    # canonical up (50x), alternate at baseline: a switch, canonical high
    @test DS.readout_switch(50.0, Union{Nothing, Float64}[BL], BL) == (true, "high", 1.0)
    # canonical not, alternate not: no switch
    sw, st, osf = DS.readout_switch(1.0, Union{Nothing, Float64}[1.01 * BL], BL)
    @test sw === false && st === nothing && osf === nothing
    # canonical down, alternate up: a switch, canonical low
    sw, st, osf = DS.readout_switch(0.2, Union{Nothing, Float64}[2.0 * BL], BL)
    @test sw === true && st == "low" && osf ≈ 2.0
    # unvalued readout: no claim
    @test DS.readout_switch(nothing, Union{Nothing, Float64}[1.0], BL) == (nothing, nothing, nothing)
    # nothing-valued overlays are skipped
    sw, _, _ = DS.readout_switch(1.0, Union{Nothing, Float64}[nothing], BL)
    @test sw === false
end

end
