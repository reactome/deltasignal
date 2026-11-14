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
    # the high-basin re-solve converges in one sweep only if already there;
    # whatever the cap leaves unconverged is counted, never served
    @test r.diagnostics["switch_sweeps_capped"] >= 0
    for c in sw_comps(r)
        @test !isempty(c["alt_activities"])
    end
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
