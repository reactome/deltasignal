# specs/022: an inhibitor that contains its own reaction's input.
#
# Fixture = the WNT shape traced in specs/021: target T has AND inputs L (the
# ligand) and F (a partner at baseline), and OR-flagged inhibitors C1, C2 whose
# containment includes L. C1 and C2 are COMPUTED from L (assembly with Y1, Y2),
# so they track it through the network, as a sequestering complex does. The
# first version of this suite pinned them to L's fold instead, and so never saw
# the case the PR #72 review found: an inhibitor that contains the input but is
# not computed from it, where the rule amplified (L = 5x read 21x).
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
node(u) = u => DS.NetworkNode(u, "R-" * u, "protein", nothing, u, BL)
edge(s, t, pos, et; and = true) = DS.LogicNetworkEdge(s, t, and, pos, 1.0, et)
const CONT = Dict("R-C1" => Set(["R-L", "R-Y1"]), "R-C2" => Set(["R-L", "R-Y2"]))

"T = L AND F, inhibited by C1 = L+Y1 and C2 = L+Y2 (assembly), optionally by C3."
function fixture(; containment = CONT, extra_inhibitor = false, tracking = true, or_inputs = false)
    nodes = Dict(node.(["T", "L", "F", "C1", "C2", "C3", "Y1", "Y2", "Z"]))
    edges = [edge("L", "T", true, "input"; and = !or_inputs),
             edge("F", "T", true, "input"; and = !or_inputs),
             edge("C1", "T", false, "regulator"; and = false),
             edge("C2", "T", false, "regulator"; and = false)]
    if tracking
        append!(edges, [edge("L", "C1", true, "assembly"), edge("Y1", "C1", true, "assembly"),
                        edge("L", "C2", true, "assembly"), edge("Y2", "C2", true, "assembly")])
    else
        # C1, C2 contain L by Reactome's containment table, but here they are
        # made from Z alone, so they do not move with L.
        append!(edges, [edge("Z", "C1", true, "assembly"), edge("Z", "C2", true, "assembly")])
    end
    extra_inhibitor && push!(edges, edge("C3", "T", false, "regulator"; and = false))
    DS.ReactionNetwork(nodes, edges, Dict{String, DS.SetExpansionMapping}(), Set{String}(), containment)
end

"Fold of T with L at `x`, F at 1, Y1 at `y1`, C3 at `c3` (UI = fold)."
function solve(net; x, y1 = 1.0, c3 = 1.0, w = "off")
    with_env("DS_SELF_INHIBITOR_WEIGHT" => w, "DS_ASSEMBLY_LIMITING" => nothing) do
        obs = Dict("L" => (x, 1.0), "F" => (1.0, 1.0), "Y1" => (y1, 1.0), "Y2" => (1.0, 1.0),
                   "Z" => (1.0, 1.0), "C3" => (c3, 1.0))
        DS.solve_steady_state(net, obs, DS.SteadyStateParams(1.0, 0.1, 500, 1e-6, "penalty"))
    end
end
t_fold(net; kw...) = (r = solve(net; kw...); r.node_activities["T"] / BL)

@testset "self-contained inhibition (specs/022)" begin

@testset "default is 0.1; off reproduces the double count exactly" begin
    with_env("DS_SELF_INHIBITOR_WEIGHT" => nothing) do
        @test DS.resolve_reaction_eval_config().self_inhibitor_weight == 0.1
    end
    with_env("DS_SELF_INHIBITOR_WEIGHT" => "off") do
        @test DS.resolve_reaction_eval_config().self_inhibitor_weight == -1.0
    end
    net = fixture()
    for (x, want) in ((0.5, 2.0), (0.2, 2.0), (0.05, 0.5), (2.0, 0.5))
        r = solve(net; x = x)
        @test r.node_activities["T"] / BL ≈ want rtol = 1e-6   # x * min(x^-2, 10)
        @test r.diagnostics["self_inhibitors"] == 0
        @test r.diagnostics["self_inhibitor_rule"] == "off"
    end
    @test t_fold(net; x = 0.5, w = nothing) ≈ 0.5^0.9 rtol = 1e-6   # unset = the default
end

@testset "with w: the input sets the direction, the inhibitors only damp it" begin
    net = fixture()
    for x in (0.05, 0.2, 0.5, 2.0, 5.0, 20.0)
        r = solve(net; x = x, w = "0.1")
        @test r.node_activities["T"] / BL ≈ x^0.9 rtol = 1e-6       # x^(1-w), however many
        @test r.diagnostics["self_inhibitors"] == 2
        @test r.diagnostics["self_inhibitor_rule"] == "on"
    end
    ups = [t_fold(net; x = x, w = "0.1") for x in (1.0, 2.0, 5.0, 20.0)]
    downs = [t_fold(net; x = x, w = "0.1") for x in (1.0, 0.5, 0.2, 0.05)]
    @test issorted(ups) && issorted(downs; rev = true)
    @test t_fold(net; x = 0.5, w = "0") ≈ 0.5 rtol = 1e-6       # w = 0: input alone
    @test t_fold(net; x = 0.5, w = "1") ≈ 1.0 rtol = 1e-6       # w = 1: the single-count cancel
    @test t_fold(net; x = 1.0, w = "0.1") ≈ 1.0 rtol = 1e-9     # baseline untouched
end

@testset "an inhibitor that contains the input but is not computed from it is left alone" begin
    # The PR #72 review's probe: with containment alone the rule divided L's
    # fold out of an inhibitor that never moved with L, and T read 5^1.9 = 21x.
    net = fixture(tracking = false)
    for x in (0.2, 5.0, 80.0)
        r = solve(net; x = x, w = "0.1")
        @test r.node_activities["T"] / BL ≈ min(x, 100.0) rtol = 1e-6   # same as off
        @test r.diagnostics["self_inhibitors"] == 0                       # not flagged
    end
    @test isempty(DS.self_contained_inhibitor_map(net))
end

@testset "the rule only ever weakens an inhibitor" begin
    net = fixture()
    # C1 doubles for its own reason (Y1 up) while L is at baseline: full strength.
    @test t_fold(net; x = 1.0, y1 = 2.0, w = "0.1") ≈ t_fold(net; x = 1.0, y1 = 2.0) rtol = 1e-9
    @test t_fold(net; x = 1.0, y1 = 2.0) ≈ 0.5 rtol = 1e-6
    # L up AND C1 up for its own reason: T may not end up further from what the
    # old model says than the damping explains -- the inhibitor is never made
    # stronger or reversed.
    for (x, y1) in ((5.0, 2.0), (0.2, 3.0), (0.2, 0.5), (5.0, 0.5))
        on, off = t_fold(net; x = x, y1 = y1, w = "0.1"), t_fold(net; x = x, y1 = y1)
        noinh = min(x, 100.0)                          # T with both inhibitors at baseline
        @test min(off, noinh) - 1e-9 <= on <= max(off, noinh) + 1e-9
    end
    # An unrelated inhibitor is untouched.
    net3 = fixture(extra_inhibitor = true)
    @test t_fold(net3; x = 1.0, c3 = 2.0, w = "0.1") ≈ 0.5 rtol = 1e-6
    @test solve(net3; x = 1.0, c3 = 2.0, w = "0.1").diagnostics["self_inhibitors"] == 2
end

@testset "a knockout at exactly 0 does not jump to full de-repression" begin
    # OR inputs, so T survives L = 0. The review found the rule then read the
    # inhibitor as 0 and de-repressed to the 10x ceiling (T = 5 against 0.5 off).
    net = fixture(or_inputs = true)
    for x in (0.0, 1e-9, 1e-6)
        on, off = t_fold(net; x = x, w = "0.1"), t_fold(net; x = x)
        @test on <= off + 1e-9                      # never MORE de-repression than off
    end
    # No jump: T moves monotonically as the knockdown deepens to exactly 0.
    # (The damped part f_s^(w/n) falls slowly -- (1e-9)^0.05 is still 0.35 -- so
    # a total loss of the shared input is not damped; with AND inputs, as every
    # catalog case has, the output is 0 anyway.)
    ts = [t_fold(net; x = x, w = "0.1") for x in (1e-1, 1e-2, 1e-3, 1e-6, 1e-9, 0.0)]
    @test issorted(ts[1:3])                        # de-repression grows as L falls ...
    @test isapprox(ts[4], ts[6]; rtol = 1e-5)      # ... and below the floor depth is moot
    @test abs(ts[3] - ts[4]) / ts[3] < 0.01        # no jump crossing the floor
    # (the 0.1% dip past the floor is L's own share of the OR input still falling
    # while the de-repression has saturated -- genuine, not a reversal)
end

@testset "inert without a containment table, and says so" begin
    r = solve(fixture(containment = Dict{String, Set{String}}()); x = 0.5, w = "0.1")
    @test r.node_activities["T"] / BL ≈ 2.0 rtol = 1e-6
    @test r.diagnostics["self_inhibitors"] == 0
    @test r.diagnostics["self_inhibitor_rule"] == "inert: no containment table"
end

@testset "invalid weights are startup errors" begin
    for bad in ("1.5", "-0.1", "abc", "NaN")
        with_env("DS_SELF_INHIBITOR_WEIGHT" => bad) do
            @test_throws ArgumentError DS.resolve_reaction_eval_config()
        end
    end
    with_env("DS_SELF_INHIBITOR_WEIGHT" => "0.1") do
        @test DS.resolve_reaction_eval_config().self_inhibitor_weight == 0.1
    end
end

@testset "containment.csv is read beside the logic network, self rows dropped" begin
    mktempdir() do d
        write(joinpath(d, "containment.csv"),
              "stable_id,contains_stable_id,reactome_release\nR-C1,R-C1,97\nR-C1,R-L,97\nR-C2,R-L,97\n")
        c = DS.parse_containment(joinpath(d, "logic_network.csv"))
        @test c == Dict("R-C1" => Set(["R-L"]), "R-C2" => Set(["R-L"]))
        @test isempty(DS.parse_containment(joinpath(d, "sub", "logic_network.csv")))
    end
end

@testset "containment survives the JSON round trip the API and CLI use" begin
    net = fixture()
    j = DS.network_containment_json(net)
    @test Set(keys(j)) == Set(["R-C1", "R-C2"])       # only stids present in the network
    @test DS.containment_from_json(j) == CONT
    @test isempty(DS.containment_from_json(nothing))
end

@testset "the map matches per target node, never a sibling's input" begin
    nodes = Dict(node.(["T1", "T2", "L1", "L2", "C"]))
    edges = [edge("L1", "T1", true, "input"), edge("L2", "T2", true, "input"),
             edge("L1", "C", true, "assembly"),
             edge("C", "T1", false, "regulator"), edge("C", "T2", false, "regulator"),
             edge("C", "L1", false, "depletion")]
    net = DS.ReactionNetwork(nodes, edges, Dict{String, DS.SetExpansionMapping}(), Set{String}(),
                             Dict("R-C" => Set(["R-L1"])))
    @test DS.self_contained_inhibitor_map(net) == Dict(("C", "T1") => ["L1"])
end

@testset "influence scores use the same model when given the network" begin
    net = fixture()
    r = solve(net; x = 0.5, w = "0.1")
    rx = DS.convert_to_reaction_network(net)
    with_env("DS_SELF_INHIBITOR_WEIGHT" => "0.1") do
        with_net = DS.compute_influence_scores(r, rx; network = net)
        without = DS.compute_influence_scores(r, rx)
        @test with_net != without        # the rule changes the local sensitivities
    end
end

end
