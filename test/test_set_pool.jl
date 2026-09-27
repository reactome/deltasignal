# specs/033: a set-valued catalyst or regulator as ONE pool node.
#
# Fixture = the RAF shape: members M1, M2, M3 of a catalyst set feed the pool P
# by `set_member` edges; P catalyses T together with the input X (T = X AND P).
# The pool's mode decides how one member's change reaches T.
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

function fixture()
    nodes = Dict(node.(["M1", "M2", "M3", "P", "X", "T", "A", "B", "O"]))
    edges = [edge("M1", "P", true, "set_member"; and = false),
             edge("M2", "P", true, "set_member"; and = false),
             edge("M3", "P", true, "set_member"; and = false),
             edge("X", "T", true, "input"), edge("P", "T", true, "catalyst"),
             # an ordinary OR node, which is NOT a pool: must not depend on the mode
             edge("A", "O", true, "output"; and = false), edge("B", "O", true, "output"; and = false)]
    DS.ReactionNetwork(nodes, edges, Dict{String, DS.SetExpansionMapping}())
end

function solve(mode; m1 = 1.0, m2 = 1.0, m3 = 1.0, a = 1.0)
    with_env("DS_SET_POOL_MODE" => mode, "DS_SELF_INHIBITOR_WEIGHT" => "off") do
        obs = Dict("M1" => (m1, 1.0), "M2" => (m2, 1.0), "M3" => (m3, 1.0),
                   "X" => (1.0, 1.0), "A" => (a, 1.0), "B" => (1.0, 1.0))
        DS.solve_steady_state(fixture(), obs, DS.SteadyStateParams(1.0, 0.1, 500, 1e-6, "penalty"))
    end
end
fold(r, u) = r.node_activities[u] / BL

@testset "set pool nodes (specs/033)" begin

@testset "set_pool_value, by hand" begin
    b = BL
    pv(folds, mode) = DS.set_pool_value([f * b for f in folds], b, mode) / b
    @test pv([2.0, 3.0, 1.0], "product") ≈ 6.0
    @test pv([0.0, 3.0, 1.0], "product") == 0.0
    @test pv([80.0, 80.0, 1.0], "product") ≈ 100.0          # capped at 100x
    @test pv([2.0, 0.25, 1.0], "extreme") ≈ 0.25            # |log .25| > |log 2|
    @test pv([2.0, 0.5, 1.0], "extreme") ≈ 0.5              # a tie goes to the lower fold
    @test pv([80.0, 1.0, 1.0], "extreme") ≈ 80.0
    @test pv([0.0, 50.0], "extreme") == 0.0                 # a knockout wins outright
    @test pv([1.0, 1.0, 1.0], "extreme") ≈ 1.0
    @test pv([8.0, 1.0, 1.0], "geomean") ≈ 2.0
    @test pv([0.0, 8.0], "geomean") == 0.0
    @test_throws ArgumentError DS.set_pool_value([b], b, "mean")
end

@testset "default is product" begin
    with_env("DS_SET_POOL_MODE" => nothing) do
        @test DS.resolve_reaction_eval_config().set_pool_mode == "product"
    end
end

@testset "one member up 2x, and one knocked out" begin
    for (mode, up, ko) in (("product", 2.0, 0.0), ("extreme", 2.0, 0.0),
                           ("geomean", 2.0^(1/3), 0.0), ("mean", 4/3, 2/3))
        r = solve(mode; m1 = 2.0)
        @test fold(r, "P") ≈ up rtol = 1e-6
        @test fold(r, "T") ≈ up rtol = 1e-6
        r = solve(mode; m1 = 0.0)
        @test fold(r, "P") ≈ ko atol = 1e-9
        @test fold(r, "T") ≈ ko atol = 1e-9
    end
end

@testset "several members up: product multiplies, extreme and geomean do not" begin
    @test fold(solve("product"; m1 = 4.0, m2 = 4.0), "P") ≈ 16.0 rtol = 1e-6
    @test fold(solve("extreme"; m1 = 4.0, m2 = 4.0), "P") ≈ 4.0 rtol = 1e-6
    @test fold(solve("geomean"; m1 = 4.0, m2 = 4.0), "P") ≈ 16.0^(1/3) rtol = 1e-6
end

@testset "an ordinary OR node does not depend on the mode" begin
    vals = [fold(solve(m; a = 3.0), "O") for m in ("product", "extreme", "geomean", "mean")]
    @test all(v -> v ≈ vals[1], vals)
    @test vals[1] ≈ 2.0 rtol = 1e-6                        # mean(3, 1): the OR rule
end

@testset "a typo in the mode is an error" begin
    for bad in ("Product", "max", "")
        with_env("DS_SET_POOL_MODE" => bad) do
            @test_throws ArgumentError DS.resolve_reaction_eval_config()
        end
    end
end

end
